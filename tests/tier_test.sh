#!/usr/bin/env bash
# tier.py on synthetic catalogues (seconds, no download): a catalogue whose sequences carry RefSeq names (a 2bit
# converted from the RefSeq copy of the assembly) must give the same tier and counts as one with GenBank names, and a
# gap on a RefSeq-named chromosome must count as a gap on a chromosome.
#
#   tests/tier_test.sh WORKDIR
set -uo pipefail
W=${1:?workdir}; HERE=$(cd "$(dirname "$0")/.." && pwd)
/bin/rm -rf "$W"; mkdir -p "$W"
fail=0
check () {  # name expected_row
    local got; got=$(python3 -I "$HERE/src/t2t_centromeres/tier.py" "$W/$1" "$W/$1/STATUS.tsv" 2>&1 | tail -1)
    if [ "$got" = "$2" ]; then echo "ok   $1"; else echo "FAIL $1: got '$got', expected '$2'"; fail=1; fi
}
make_case () {  # name naming(genbank|refseq) gap(yes|no)
    local D="$W/$1" a=GCA_000000001.1
    mkdir -p "$D/inputs/$a" "$D/satellome/$a/reports"
    printf 'accession\tstatus\n%s\tcalled\n' "$a" > "$D/STATUS.tsv"
    {
        echo "# Assembly name:  test"
        printf '# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n'
        printf 'chr1\tassembled-molecule\t1\tChromosome\tCM000001.1\t=\tNC_000001.1\tPrimary Assembly\t1000000\tna\n'
        printf 'chr2\tassembled-molecule\t2\tChromosome\tCM000002.1\t=\tNC_000002.1\tPrimary Assembly\t900000\tna\n'
        printf 'MT\tassembled-molecule\tMT\tMitochondrion\tCM000003.1\t=\tNC_000003.1\tnon-nuclear\t20000\tna\n'
    } > "$D/inputs/$a/$a.assembly_report.txt"
    if [ "$2" = refseq ]; then c1=NC_000001.1; c2=NC_000002.1; else c1=CM000001.1; c2=CM000002.1; fi
    {
        echo "# telomere-check v0.1.0"
        printf 'chr\tlength\tleft_count\tleft_density\tleft_motif\tleft_status\tright_count\tright_density\tright_motif\tright_status\n'
        printf '%s\t1000000\t1600\t0.96\tCCCTAA\tPRESENT\t1600\t0.96\tTTAGGG\tPRESENT\n' "$c1"
        printf '%s\t900000\t1600\t0.96\tCCCTAA\tPRESENT\t1600\t0.96\tTTAGGG\tPRESENT\n' "$c2"
    } > "$D/satellome/$a/reports/telomeres.tsv"
    : > "$D/satellome/$a/GCA_000000001v1.gaps.bed"
    [ "$3" = yes ] && printf '%s\t5000\t5100\n' "$c2" > "$D/satellome/$a/GCA_000000001v1.gaps.bed"
}
T=$'\t'
make_case genbank_t2t genbank no
check genbank_t2t "GCA_000000001.1${T}2${T}0${T}0${T}0${T}2${T}0${T}0${T}real_t2t"
make_case refseq_t2t refseq no
check refseq_t2t "GCA_000000001.1${T}2${T}0${T}0${T}0${T}2${T}0${T}0${T}real_t2t"
make_case genbank_gap genbank yes
check genbank_gap "GCA_000000001.1${T}2${T}1${T}1${T}0${T}2${T}0${T}0${T}gapped"
make_case refseq_gap refseq yes
check refseq_gap "GCA_000000001.1${T}2${T}1${T}1${T}0${T}2${T}0${T}0${T}gapped"
# a chromosome with no telomere row under either name is still refused
make_case refseq_missing refseq no
sed -i.bak '/^NC_000002.1/d' "$W/refseq_missing/satellome/GCA_000000001.1/reports/telomeres.tsv"
if python3 -I "$HERE/src/t2t_centromeres/tier.py" "$W/refseq_missing" "$W/refseq_missing/STATUS.tsv" > /dev/null 2> "$W/refseq_missing.err" \
    || ! grep -q 'CM000002.1 has no telomere row' "$W/refseq_missing.err"; then
    echo "FAIL refseq_missing: not refused with the chromosome named"; fail=1
else echo "ok   refseq_missing"; fi
[ "$fail" -eq 0 ] && echo "tier test: all passed"
exit "$fail"
