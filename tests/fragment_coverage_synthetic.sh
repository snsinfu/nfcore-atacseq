#!/usr/bin/env bash
#
# Regression test for paired-end coverage in modules/local/bedtools_genomecov.nf.
#
# `bedtools genomecov -pc` silently drops any pair whose reverse-strand mate
# starts upstream of its forward-strand mate (dovetailing pairs, i.e. fragment
# length < read length). This builds a tiny coordinate-sorted paired-end BAM of
# such pairs and checks that:
#   1. `bedtools genomecov -pc` covers 0 bases (the regression), and
#   2. building one fragment per pair from its positive-TLEN mate covers the
#      expected number of bases.
#
# Requires: samtools, bedtools, awk, sort.
# Run: bash tests/fragment_coverage_synthetic.sh

set -euo pipefail

READ_LEN=125
INSERT=60
N=2000
SPACING=200
CONTIG=1000000
CHR=chrBig
READ1_START=1000

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# --- build synthetic dovetailing PE BAM -------------------------------------
awk -v rl="$READ_LEN" -v ins="$INSERT" -v n="$N" -v sp="$SPACING" \
    -v chr="$CHR" -v c="$CONTIG" -v s0="$READ1_START" '
BEGIN {
    q = ""; a = ""
    for (j = 0; j < rl; j++) { q = q "I"; a = a "A" }
    print "@HD\tVN:1.6\tSO:coordinate"
    print "@SQ\tSN:" chr "\tLN:" c
    for (i = 0; i < n; i++) {
        p = s0 + i * sp            # forward read start
        m = p + ins - rl           # reverse read start (upstream: m < p)
        printf("f%d\t99\t%s\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s\n",  i, chr, p, rl, m,  ins, a, q)
        printf("f%d\t147\t%s\t%d\t60\t%dM\t=\t%d\t%d\t%s\t%s\n", i, chr, m, rl, p, -ins, a, q)
    }
}' > "$tmp/dovetail.sam"

samtools view -b "$tmp/dovetail.sam" | samtools sort -o "$tmp/dovetail.bam" -
printf '%s\t%d\n' "$CHR" "$CONTIG" > "$tmp/genome.sizes"

covered() {  # sum covered bases from a bedGraph stream
    awk '{s += $3 - $2} END {print s + 0}'
}

# --- old recipe: bedtools genomecov -pc -------------------------------------
old="$(bedtools genomecov -ibam "$tmp/dovetail.bam" -bg -pc | covered)"

# --- new recipe: positive-TLEN fragment BED ---------------------------------
new="$(samtools view -f 0x2 -F 0x904 "$tmp/dovetail.bam" \
    | awk 'BEGIN{OFS="\t"} $9 > 0 {print $3, $4 - 1, $4 - 1 + $9}' \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bedtools genomecov -i - -g "$tmp/genome.sizes" -bg | covered)"

expected=$(( N * INSERT ))

echo "bedtools genomecov -pc covered bases : $old   (expected 0; -pc drops dovetailing pairs)"
echo "one fragment per pair covered bases  : $new   (expected $expected)"

fail=0
if [ "$old" -ne 0 ]; then
    echo "FAIL: bedtools genomecov -pc covered $old bases; expected 0" >&2
    fail=1
fi
if [ "$new" -ne "$expected" ]; then
    echo "FAIL: fragment coverage was $new bases; expected $expected" >&2
    fail=1
fi
[ "$fail" -eq 0 ] && echo "PASS"
exit "$fail"
