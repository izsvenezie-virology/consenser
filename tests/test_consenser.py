import io

from consenser.consenser import (
    Mutation,
    Thresholds,
    apply_indel,
    fasta_format,
    get_indels,
    get_snp,
    read_fasta,
    read_vcf,
)


def test_read_fasta_multiple_sequences():
    content = ">chr1\nACGT\n>chr2\nTTAA"
    reference = read_fasta(io.StringIO(content))

    assert reference == {"chr1": "ACGT", "chr2": "TTAA"}


def test_read_vcf_parses_mutations():
    content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
PB2	33	.	A	G	49314	PASS	DP=2201;AF=1.000000;SB=0;DP4=0,0,1864,337
"""
    vcf = read_vcf(io.StringIO(content))

    muts = vcf["PB2"]["33"]
    assert len(muts) == 1
    mut = muts[0]
    assert mut == Mutation("PB2", 33, "A", "G", 1.0)


def test_fasta_format_wraps_width():
    seq = "A" * 10
    assert fasta_format(seq) == [seq]


def test_get_snp_returns_iupac_for_supported_variants():
    muts = [
        Mutation("chr1", 1, "A", "G", 0.3),
        Mutation("chr1", 1, "A", "C", 0.3),
    ]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "V"


def test_get_snp_uses_highest_when_below_threshold():
    muts = [
        Mutation("chr1", 1, "A", "G", 0.1),
        Mutation("chr1", 1, "A", "C", 0.1),
    ]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "A"


def test_apply_indel_insertion():
    seq = list("ACG")
    indel = Mutation("chr1", 2, "C", "CT", 0.8)

    assert apply_indel(seq, indel) == ["A", "CT", "G"]


def test_apply_indel_deletion():
    seq = list("ACGT")
    indel = Mutation("chr1", 2, "CG", "C", 0.8)

    assert apply_indel(seq, indel) == ["A", "C", "", "T"]


def test_get_indels_empty_returns_none():
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_indels([], thresholds) is None
