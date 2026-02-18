import io
from decimal import Decimal

import pytest

from consenser.consenser import (
    MalformedFasta,
    Mutation,
    Thresholds,
    apply_indel,
    create_consensus_sequence,
    fasta_format,
    get_indels,
    get_snp,
    read_coverage,
    read_fasta,
    read_vcf,
)


def test_read_fasta_single_sequence():
    content = ">chr1\nACGT"
    reference = read_fasta(io.StringIO(content))

    assert reference == {"chr1": "ACGT"}


def test_read_fasta_multiple_sequences():
    content = ">chr1\nACGT\n>chr2\nTTAA"
    reference = read_fasta(io.StringIO(content))

    assert reference == {"chr1": "ACGT", "chr2": "TTAA"}


def test_read_fasta_malformed():
    content = "ACTG\n>chr1\nACGT"
    with pytest.raises(MalformedFasta):
        read_fasta(io.StringIO(content))


def test_read_vcf_parses_mutation():
    content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
PB2	33	.	A	G	49314	PASS	DP=2201;AF=1.000000;SB=0;DP4=0,0,1864,337
"""
    vcf = read_vcf(io.StringIO(content))

    muts = vcf["PB2"][33]
    assert len(muts) == 1
    mut = muts[0]
    assert mut == Mutation("PB2", 33, "A", "G", 1.0)


def test_read_vcf_parses_mutations():
    content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
PB2	33	.	A	G	49314	PASS	DP=2201;AF=0.500000;SB=0;DP4=0,0,1864,337
PB2	33	.	A	T	49314	PASS	DP=2201;AF=0.499999;SB=0;DP4=0,0,1864,337
PB2	33	.	A	C	49314	PASS	DP=2201;AF=0.000001;SB=0;DP4=0,0,1864,337
"""
    vcf = read_vcf(io.StringIO(content))

    muts = vcf["PB2"][33]
    assert len(muts) == 3
    assert muts[0] == Mutation("PB2", 33, "A", "G", Decimal("0.5"))
    assert muts[1] == Mutation("PB2", 33, "A", "T", Decimal("0.499999"))
    assert muts[2] == Mutation("PB2", 33, "A", "C", Decimal("0.000001"))


def test_read_coverage_file():
    content = """PB2	1	4
PB2	2	5
PB2	3	6
PB2	4	8
PB2	5	8
PB2	6	8
"""
    coverage = read_coverage(io.StringIO(content))

    assert coverage == {"PB2": [4, 5, 6, 8, 8, 8]}


def test_get_snp_over_max_threshold():
    muts = [Mutation("chr1", 1, "A", "G", Decimal("0.76"))]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "G"


def test_get_snp_equal_max_threshold():
    muts = [Mutation("chr1", 1, "A", "G", Decimal("0.75"))]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "R"


def test_get_snp_below_max_threshold():
    muts = [Mutation("chr1", 1, "A", "G", Decimal("0.74"))]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "R"


def test_get_snp_over_min_threshold():
    muts = [Mutation("chr1", 1, "A", "G", Decimal("0.26"))]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "R"


def test_get_snp_equal_min_threshold():
    muts = [Mutation("chr1", 1, "A", "G", Decimal("0.25"))]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "R"


def test_get_snp_below_min_threshold():
    muts = [Mutation("chr1", 1, "A", "G", Decimal("0.24"))]
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "A"


def test_get_most_present_snp_ref():
    muts = [
        Mutation("chr1", 1, "A", "G", Decimal("0.33")),
        Mutation("chr1", 1, "A", "T", Decimal("0.33")),
    ]
    thresholds = Thresholds(lower=1.0, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "A"


def test_get_most_present_snp_alt():
    muts = [
        Mutation("chr1", 1, "A", "G", Decimal("0.33")),
        Mutation("chr1", 1, "A", "T", Decimal("0.34")),
    ]
    thresholds = Thresholds(lower=1.0, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "T"


def test_get_most_present_snp_fifty_fifty():
    muts = [
        Mutation("chr1", 1, "A", "G", Decimal("0.5")),
    ]
    thresholds = Thresholds(lower=1.0, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "R"


def test_get_most_present_snp_thirds():
    muts = [
        Mutation("chr1", 1, "A", "G", Decimal("0.33")),
        Mutation("chr1", 1, "A", "T", Decimal("0.33")),
        Mutation("chr1", 1, "A", "C", Decimal("0.01")),
    ]
    thresholds = Thresholds(lower=1.0, indel=0.5, coverage=10)

    assert get_snp(muts, thresholds) == "D"


def test_get_most_frequent_indel():
    muts = [
        Mutation("chr1", 1, "A", "ATCG", Decimal("0.33")),
        Mutation("chr1", 1, "A", "ATCG", Decimal("0.33")),
        Mutation("chr1", 1, "A", "ATC", Decimal("0.34")),
    ]
    thresholds = Thresholds(lower=0.25, indel=0.25, coverage=10)

    assert get_indels(muts, thresholds) == muts[2]


def test_no_indels_thresholds():
    muts = [
        Mutation("chr1", 1, "A", "ATCG", Decimal("0.33")),
        Mutation("chr1", 1, "A", "ATCG", Decimal("0.33")),
        Mutation("chr1", 1, "A", "ATC", Decimal("0.34")),
    ]
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)

    assert get_indels(muts, thresholds) is None


def test_get_indels_empty():
    thresholds = Thresholds(lower=0.25, indel=0.5, coverage=10)

    assert get_indels([], thresholds) is None


def test_indels_snp():
    muts = [
        Mutation("chr1", 1, "A", "ATCG", Decimal("0.8")),
        Mutation("chr1", 1, "A", "C", Decimal("0.6")),
        Mutation("chr1", 1, "A", "T", Decimal("0.1")),
    ]
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)

    assert get_indels(muts, thresholds) == muts[0]
    assert get_snp(muts, thresholds) == "M"


def test_apply_indel_insertion():
    seq = list("ACG")
    indel = Mutation("chr1", 2, "C", "CT", 0.8)

    assert apply_indel(seq, indel) == ["A", "CT", "G"]


def test_apply_indel_deletion():
    seq = list("ACGT")
    indel = Mutation("chr1", 2, "CG", "C", 0.8)

    assert apply_indel(seq, indel) == ["A", "C", "", "T"]


def test_create_consensus_complete():
    seq = "ACTGACTGACTG"
    muts = {
        1: [
            Mutation("chr1", 1, "A", "ATCG", Decimal("0.8")),
            Mutation("chr1", 1, "A", "C", Decimal("0.6")),
            Mutation("chr1", 1, "A", "T", Decimal("0.1")),
        ],
        5: [
            Mutation("chr1", 5, "ACTG", "A", Decimal("0.8")),
            Mutation("chr1", 5, "A", "T", Decimal("0.6")),
            Mutation("chr1", 5, "A", "C", Decimal("0.1")),
        ],
    }
    cov = [10, 11, 11, 10, 11, 11, 10, 11, 11, 10, 9, 8]
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)
    consensus = create_consensus_sequence(seq, muts, cov, thresholds)
    assert consensus == "MTCGCTGWACNN"


def test_create_consensus_no_coverage():
    seq = "ACTGACTGACTG"
    muts = {
        1: [
            Mutation("chr1", 1, "A", "ATCG", Decimal("0.8")),
            Mutation("chr1", 1, "A", "C", Decimal("0.6")),
            Mutation("chr1", 1, "A", "T", Decimal("0.1")),
        ]
    }
    cov = []
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)
    consensus = create_consensus_sequence(seq, muts, cov, thresholds)
    assert consensus == "MTCGCTGACTGACTG"


def test_create_consensus_del_low_coverage():
    seq = "ACTGACTGACTG"
    muts = {1: [Mutation("chr1", 1, "ACTG", "A", Decimal("0.8"))]}
    cov = [9, 11, 11, 10, 11, 11, 10, 11, 11, 10, 9, 8]
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)
    consensus = create_consensus_sequence(seq, muts, cov, thresholds)
    assert consensus == "NCTGACTGACNN"


def test_create_consensus_ins_low_coverage():
    seq = "ACTGACTGACTG"
    muts = {1: [Mutation("chr1", 1, "A", "ATCG", Decimal("0.8"))]}
    cov = [9, 11, 11, 10, 11, 11, 10, 11, 11, 10, 9, 8]
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)
    consensus = create_consensus_sequence(seq, muts, cov, thresholds)
    assert consensus == "NCTGACTGACNN"


def test_create_consensus_snp_low_coverage():
    seq = "ACTGACTGACTG"
    muts = {1: [Mutation("chr1", 1, "A", "T", Decimal("0.8"))]}
    cov = [9, 11, 11, 10, 11, 11, 10, 11, 11, 10, 9, 8]
    thresholds = Thresholds(lower=0.25, indel=0.50, coverage=10)
    consensus = create_consensus_sequence(seq, muts, cov, thresholds)
    assert consensus == "NCTGACTGACNN"


def test_fasta_format():
    sequence = "01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789"
    fasta = fasta_format(sequence)

    assert len(fasta) == 3
    assert (
        fasta[0]
        == fasta[1]
        == "0123456789012345678901234567890123456789012345678901234567890123456789"
    )
    assert fasta[2] == "012345678901234567890123456789"
