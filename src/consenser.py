import click
from click.exceptions import ClickException
from click.types import File, Path
from copy import deepcopy


class Mutation():
    def __init__(self, chrom: str, pos: int, ref: str, alt: str, freq: float) -> None:
        self.chromosome: str = chrom
        self.position: int = pos
        self.reference: str = ref
        self.alteration: str = alt
        self.frequency: float = freq

    def is_indel(self) -> bool:
        return len(self.reference) > 1 or len(self.alteration) > 1

    def is_deletion(self) -> bool:
        return len(self.reference) > 1

    def is_insertion(self) -> bool:
        return len(self.alteration) > 1

    def __str__(self) -> str:
        return(f'Pos: {self.position}; Ref: {self.reference}; Alt: {self.alteration}; Freq: {self.frequency}')


@click.command()
@click.option('-c', '--coverage', type=File('r'), required=True, help='The coverage file, a tab separated file with Chrom, Position and Coverage columns, without header.')
@click.option('-r', '--reference', type=File('r'), required=True, help='The reference file in fasta format.')
@click.option('-v', '--vcf', type=File('r'), required=True, help='The VCF file. This file must contain the allele frequency ("AF=") in the INFO column.')
@click.option('-o', '--output', type=File('w'), default='-', help='The output file. [default: stout]')
@click.option('-d', '--deg', type=float, nargs=2, help='The upper and lower limit to insert a degeneration; in percentage.')
@click.option('-n', '--no-deg', type=float, nargs=1, help='The minimum AF percentage to consider a snp. It avoids degenerations.')
@click.option('-m', '--min-cov', type=int, default=10, show_default=True, help='Minimum coverage to not mask a base.')
@click.option('-w', '--width', type=int, default=70, show_default=True, help='The width of the Fasta files in output.')
@click.option('-s', '--split', type=Path(), help='Creates a file for each one of the sequences. The keyword "chromHere" in the path will be replaced with the original sequence name.')
@click.option('-a', '--alter_names', type=str, help='Replace the sequence name. The keyword "chrmoHere" will be replaced with the original sequence name.')
@click.option('--indels-lim', type=float, default=50.0, show_default=True, help='Set the minimum limit to consider an indel.')
def cli(coverage: File, reference: File, vcf: File, output: File, deg: tuple[float, float], no_deg: float,
        min_cov: int, width: int, split: Path, alter_names: str, indels_lim: float):
    '''Creates a consensus sequence from a reference and the VCF file. Low coverage regions are masked using the coverage file.'''
    af_lims = parse_limits(deg, no_deg, indels_lim)
    cov_by_chrom = read_coverage(coverage, min_cov)
    ref_by_chrom = read_reference(reference)
    vcf_by_chrom = read_vcf(vcf)
    chroms = ref_by_chrom.keys()

    consensus = {}

    for chrom in chroms:
        cons_seq = create_consensus_sequence(
            ref_by_chrom[chrom], vcf_by_chrom[chrom], cov_by_chrom[chrom], af_lims)
        cons_name = chrom
        if alter_names:
            cons_name = alter_names.replace("chromHere", chrom)
        consensus[chrom] = (cons_name, cons_seq)

    write_consensus(consensus, output, split, width)


def parse_limits(deg: tuple[float, float], no_deg: float, indels_lim: float) -> tuple[float, float, float]:
    '''Perfrorms checks and setup limits'''
    if not deg and not no_deg:
        raise ClickException("'-d'/'--deg' or '-n'/'--no-deg' is required")
    if deg and no_deg:
        raise ClickException("Use '-d'/'--deg' or '-n'/'--no-deg', not both")
    if (no_deg):
        deg = (no_deg, no_deg)
    if deg[1] < 50:
        raise ClickException("The upper limit must be at least 50%")
    return deg + (indels_lim,)


def read_coverage(cov_file: File, min_cov: int) -> dict[str, list[int]]:
    '''Reads the coverage file and returns low coverage positions grouped by chromosome'''
    chroms_cov: dict[str, list[int]] = {}
    for line in cov_file.read().split('\n'):
        if line == '':
            continue
        chrom, pos, cov = line.split('\t')
        if chrom not in chroms_cov:
            chroms_cov[chrom] = []
        if int(cov) < min_cov:
            chroms_cov[chrom].append(int(pos))
    return chroms_cov


def read_reference(ref_file: File) -> dict[str, str]:
    chroms_ref = {}
    for line in ref_file.read().split('\n'):
        if line == '':
            continue
        if line.startswith('>'):
            chrom_name = line[1:]
            chroms_ref[chrom_name] = ''
            continue
        chroms_ref[chrom_name] += line
    return chroms_ref


def read_vcf(vcf_file: File) -> tuple[dict[str, dict[int, list[Mutation]]]]:
    chroms_vcf: dict[str, dict[int, list[Mutation]]] = {}
    for line in vcf_file.read().split('\n'):
        if line == '' or line.startswith('#'):
            continue
        chrom, pos, _, ref, alt, _, _, info = line.split('\t')
        pos = int(pos)
        infos = info.split(';')
        af = [i for i in infos if i.startswith('AF=')][0][3:]
        af = float(af)*100
        mut = Mutation(chrom, pos, ref, alt, af)

        if chrom not in chroms_vcf:
            chroms_vcf[chrom] = {}
        if pos not in chroms_vcf[chrom]:
            chroms_vcf[chrom][pos] = []

        chroms_vcf[chrom][pos].append(mut)
    return chroms_vcf


def create_consensus_sequence(chrom_seq: str, vcf: dict[int, list[Mutation]], cov: list, lims: tuple[float, float]):
    seq = list(chrom_seq)

    for index in range(len(seq)):
        pos = index + 1
        if seq[index] == '#':
            continue
        if pos in vcf:
            muts = vcf[pos]
            check_reference(seq, muts[0])
            snp, indel = parse_muts(muts, lims)
            if snp:
                seq[index] = snp.alteration
            if indel:
                seq = apply_indel(seq, indel)
        if pos in cov:
            seq[index] = 'N'
    return ''.join(seq).replace('#', '')


def parse_muts(muts: list[Mutation], lims: tuple[float, float, float]) -> tuple[Mutation, Mutation]:
    snp = get_consensus_snp(muts, lims[:-1])
    indel = get_consensus_indel(muts, lims[-1])
    return snp, indel


def check_reference(seq: list[str], mut: Mutation) -> None:
    '''Checks if reference and mutation's reference match'''
    index = mut.position - 1
    if seq[index] != mut.reference[0]:
        msg = f'Mismatch between mutation and reference at position {mut.position}'
        raise ValueError(msg)


def get_consensus_indel(muts: list[Mutation], freq_min: float) -> Mutation:
    '''Returns the indel to insert/delete from consensus sequence'''
    indels = [m for m in muts if m.is_indel() and m.frequency >= freq_min]
    indels.sort(key=lambda m: m.frequency, reverse=True)
    return indels[0] if indels else None


def get_consensus_snp(muts: list[Mutation], lims: tuple[float, float]) -> Mutation:
    '''Creates the snp that best fits parameters'''
    min_freq, max_freq = lims
    snps = [m for m in muts if not m.is_indel()]
    degs: list[Mutation] = []

    if not snps:
        return None

    # Create a Mutation class to use with reference nucleotide
    ref = deepcopy(snps[0])
    ref.frequency = 100.0
    ref.alteration = ref.reference

    for snp in snps:
        if snps[0].frequency >= max_freq:
            return snp
        if snp.frequency >= min_freq:
            degs.append(snp)
        ref.frequency -= snp.frequency

    if ref.frequency >= max_freq:
        return None
    if ref.frequency >= min_freq:
        degs.append(ref)

    return degenerate(degs)


def degenerate(degs: list[Mutation]) -> Mutation:
    m = degs[0]
    result = Mutation(m.chromosome, m.position, m.reference, '', 0)

    for mut in degs:
        result.alteration += mut.alteration
        result.frequency += mut.frequency
    result.alteration = to_iupac(result.alteration)
    return result


def to_iupac(nucls: str) -> str:
    sorted_nucls = ''.join(sorted(nucls))
    return iupac_names[sorted_nucls]


def apply_indel(seq: list[str], indel: Mutation) -> list[str]:
    if indel.is_insertion():
        index = indel.position - 1
        seq[index] += indel.alteration[1:]
        return seq
    for i in range(len(indel.reference) - 1):
        seq[indel.position + i] = '#'
    return seq


iupac_names = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y',  'GT': 'K',
    'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B',
    'ACGT': 'N'}


def write_consensus(cons: dict[str, (str, str)], out: File, split: str, width: int):
    for chrom in cons:
        name, seq = cons[chrom]
        name = f'>{name}'
        fasta = fasta_format(seq, width)
        fasta.insert(0, name)

        write_fasta(fasta, out)
        if split:
            f_path = split.replace('chromHere', chrom)
            with open(f_path, 'w') as f:
                write_fasta(fasta, f)


def write_fasta(lines: list[str], out:File) -> None:
    for line in lines:
        out.write(f'{line}\n')


def fasta_format(seq: str, width: int) -> list[str]:
    return [seq[i:width+i] for i in range(0, len(seq), width)]


if __name__ == '__main__':
    cli()
