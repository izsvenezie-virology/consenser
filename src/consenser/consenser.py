#! /usr/bin/env python3

from collections import defaultdict, namedtuple
from decimal import Decimal
from typing import Dict, List, Optional, TextIO

import click
from click.types import File

# fmt: off
WIDTH = 70
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "M", "AG": "R", "AT": "W",
    "CG": "S", "CT": "Y", "GT": "K",
    "ACG": "V", "ACT": "H", 
    "AGT": "D", "CGT": "B",
    "ACGT": "N",
}
# fmt: on

__author__ = "EdoardoGiussani"
__contact__ = "egiussani@izsvenezie.it"
__version__ = "1.0.0"

Thresholds = namedtuple("Thresholds", ["lower", "indel", "coverage"])
Mutation = namedtuple("Mutation", ["chrom", "pos", "ref", "alt", "freq"])


class MalformedFasta(Exception):
    def __init__(self) -> None:
        message = "Reference Fasta file is not correctly formatted"
        super().__init__(message)


def main(
    reference_file: TextIO,
    vcf_file: TextIO,
    output: TextIO,
    coverage_file: TextIO,
    header: str,
    minimum_coverage: int,
    snp_threshold: float,
    indel_threshold: float,
):
    sequences = read_fasta(reference_file)
    variants = read_vcf(vcf_file)
    coverage = read_coverage(coverage_file)

    thresholds = Thresholds(snp_threshold, indel_threshold, minimum_coverage)

    consensus: Dict[str, str] = {}

    for chrom in sequences.keys():
        cons_seq = create_consensus_sequence(
            sequences.get(chrom, ""),
            variants.get(chrom, {}),
            coverage.get(chrom, []),
            thresholds,
        )
        cons_name = header.replace("CHROM", chrom)
        consensus[cons_name] = cons_seq

    write_consensus(consensus, output)


def read_coverage(cov_file: TextIO) -> Dict[str, List[int]]:
    chroms_cov: Dict[str, List[int]] = defaultdict(list)
    if not cov_file:
        return chroms_cov
    for line in cov_file.read().split("\n"):
        if not line:
            continue
        chrom, pos, cov = line.split("\t")
        if "e" in cov:
            cov = 1_000_000_000
        chroms_cov[chrom].append(int(cov))
    return chroms_cov


def read_fasta(reference_file: TextIO) -> Dict[str, str]:
    reference = {}
    for line in reference_file.read().split("\n"):
        if not line:
            continue
        if line.startswith(">"):
            chrom_name = line[1:]
            reference[chrom_name] = ""
            continue
        try:
            reference[chrom_name] += line  # type: ignore
        except UnboundLocalError:
            raise MalformedFasta() from None
    return reference


def read_vcf(vcf_file: TextIO) -> Dict[str, Dict[int, List[Mutation]]]:
    vcf: Dict[str, Dict[int, List[Mutation]]] = defaultdict(lambda: defaultdict(list))
    for line in vcf_file.read().split("\n"):
        if not line or line.startswith("#"):
            continue

        chrom, pos, _, ref, alt, _, _, info = line.split("\t")
        pos = int(pos)
        af = [i for i in info.split(";") if i.startswith("AF=")][0][3:]

        mut = Mutation(chrom, pos, ref, alt, Decimal(af))
        vcf[chrom][pos].append(mut)
    return vcf


def create_consensus_sequence(
    chrom_seq: str,
    vcf: Dict[int, List[Mutation]],
    cov: List[int],
    thresholds: Thresholds,
):
    seq = list(chrom_seq)

    for index in range(len(seq) - 1, -1, -1):
        if cov and cov[index] < thresholds.coverage:
            seq[index] = "N"
            continue

        variants = vcf.get(index + 1, [])
        snp = get_snp(variants, thresholds)
        if snp:
            seq[index] = snp

        indel = get_indels(variants, thresholds)
        if indel:
            seq = apply_indel(seq, indel)

    return "".join(seq)


def get_indels(muts: List[Mutation], thresholds: Thresholds) -> Optional[Mutation]:
    """Returns the indel to insert/delete from consensus sequence"""
    indels = [m for m in muts if len(m.ref + m.alt) > 2 and m.freq >= thresholds.indel]
    indels.sort(key=lambda m: m.freq, reverse=True)
    return indels[0] if indels else None


def get_snp(muts: List[Mutation], thresholds: Thresholds) -> Optional[str]:
    """Creates the snp that best fits parameters"""
    snps = [m for m in muts if len(m.ref + m.alt) == 2]
    if not snps:
        return None

    nucleotides = []
    ref_freq = 1
    higher_freq = 0
    higher_nucl = ""

    for snp in snps:
        ref_freq -= snp.freq
        if snp.freq >= thresholds.lower:
            nucleotides.append(snp.alt)
        if snp.freq == higher_freq:
            higher_nucl += snp.alt
        if snp.freq > higher_freq:
            higher_nucl = snp.alt
            higher_freq = snp.freq

    if ref_freq >= thresholds.lower:
        nucleotides.append(snps[0].ref)
    if ref_freq == higher_freq:
        higher_nucl += snps[0].ref
    if ref_freq > higher_freq:
        higher_nucl = snps[0].ref

    if not nucleotides:
        nucleotides = higher_nucl

    return IUPAC["".join(sorted(nucleotides))]


def apply_indel(seq: List[str], indel: Mutation) -> List[str]:
    """Insert or delete the indel in the consensus"""
    if len(indel.alt) > 1:
        seq[indel.pos - 1] += indel.alt[1:]
        return seq
    for i in range(len(indel.ref) - 1):
        seq[indel.pos + i] = ""
    return seq


def write_consensus(consensus: Dict[str, str], out: TextIO) -> None:
    """Write the consensus sequences to Fasta file"""
    for header, sequence in consensus.items():
        out.write(f">{header}\n")
        for line in fasta_format(sequence):
            out.write(f"{line}\n")


def fasta_format(seq: str) -> List[str]:
    """Get a sequence and returns a list of subsequence of specified len"""
    return [seq[i : i + WIDTH] for i in range(0, len(seq), WIDTH)]


# fmt: off
@click.command()
@click.help_option("-h", "--help")
@click.version_option(__version__, "-v", "--version", message=f"%(prog)s, version %(version)s, by {__author__} ({__contact__})")
@click.option("-o", "--output", type=File("w"), default="-", help="The output file. [default: stout]")
@click.option("-c", "--coverage-file", type=File("r"), default=None, help="The coverage file, a tab separated file with Chrom, Position and Coverage columns, without header.")
@click.option("-H", "--header", type=str, default="CHROM", show_default=True, help='Replace the sequence name. The keyword "CHROM" will be replaced with the original sequence name.')
@click.option("-m", "--minimum-coverage", type=int, default=10, show_default=True, help="Minimum coverage to not mask a base.")
@click.option("-t", "--snp-threshold", type=Decimal, default=0.25, show_default=True, help="Minimum SNP frequency to be considered.")
@click.option("-i", "--indel-threshold", type=Decimal, default=0.50, show_default=True, help="Minimum INDELs frequency to be considered.")
@click.argument("reference_file", type=File("r"))
@click.argument("vcf_file", type=File("r"))
def cli(**kwargs):
    """Creates a consensus sequence from the reference and the VCF file.
    \b
    REFERENCE           Reference used during alignment in Fasta format.\b
    VCF                 VCF file."""
    main(**kwargs)

if __name__ == "__main__":
    cli()
