import click
from click.exceptions import ClickException
from click.types import File, Path


@click.command()
@click.option('-c', '--coverage', type=File('r'), required=True)
@click.option('-r', '--reference', type=File('r'), required=True)
@click.option('-v', '--vcf', type=File('r'), required=True)
@click.option('-d', '--deg', type=float, nargs=2)
@click.option('-n', '--no-deg', type=float, nargs=1)
@click.option('-m', '--min-cov', type=int, default=10)
@click.option('-w', '--width', type=int, default=70)
@click.option('-s', '--split', type=Path())
@click.option('-a', '--alter_names', type=str)
@click.option('--indels-lim', type=float, default=50.0)
def cli(coverage: File, reference: File, vcf: File, deg: tuple[float, float], no_deg: float,
        min_cov: int, width: int, split: Path, alter_names: str, indels_lim: float):
    if not deg and not no_deg:
        raise ClickException("'-d'/'--deg' or '-n'/--no-deg' is required")
    if deg and no_deg:
        raise ClickException("Use '-d'/'--deg' or '-n'/--no-deg', not both")
    if (no_deg):
        deg = (no_deg, no_deg)
    af_lims = deg + (indels_lim,)

    cov_by_chrom = read_coverage(coverage, min_cov)
    ref_by_chrom = read_reference(reference)
    vcf_by_chrom = read_vcf(vcf)
    chroms = ref_by_chrom.keys()

    for chrom in chroms:
        cons_seq = create_consensus_sequence(
            ref_by_chrom[chrom], vcf_by_chrom[chrom], cov_by_chrom[chrom], af_lims)


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
    chroms_snps: dict[str, dict[int, list[Mutation]]] = {}
    chroms_indels: dict[str, dict[int, list[Mutation]]] = {}
    for line in vcf_file.read().split('\n'):
        if line == '' or line.startswith('#'):
            continue
        chrom, pos, _, ref, alt, _, _, info = line.split('\t')
        infos = info.split(';')
        af = [i for i in infos if i.startswith('AF=')][0][3:]
        af = float(af)*100
        mut = Mutation(chrom, int(pos), ref, alt, af)

        mut_dict = chroms_snps
        if mut.is_indel():
            mut_dict = chroms_indels

        if chrom not in chroms_snps:
            mut_dict[chrom] = {}
        if pos not in chroms_snps[chrom]:
            mut_dict[chrom][pos] = []

        mut_dict[chrom][pos].append(mut)
    return chroms_snps, chroms_indels


def create_consensus_sequence(seq: str, vcf: dict[int, list[Mutation]], cov: list, lims: tuple[float, float]):
    mutable_seq = list(seq)

    snps = select_snps(vcf, lims[:-1])
    snps_seq = apply_snps(mutable_seq, snps)

    cov_seq = mask_low_cov(snps_seq, cov)

    indels = select_indels(vcf, lims[-1])
    indels_seq = apply_indel(cov_seq, indels)
    return ''.join(indels_seq)


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


def get_consensus_snps(muts: list[Mutation], freq_min: float) -> list[Mutation]:
    '''Returns the mutations over the minimum frequency'''
    return [m for m in muts if not m.is_indel() and m.frequency >= freq_min]


def mask_low_cov(seq: list[str], cov: list[int]):
    for pos in cov:
        seq[pos - 1] = 'N'
    return seq


def apply_snps(seq: list[str], snps: list[Mutation]) -> list[str]:
    for mut in snps:
        check_reference(seq, mut)
        seq[mut.position - 1] = to_iupac(mut.alteration)
    return seq


def apply_indel(seq: list[str], indel: Mutation) -> list[str]:
    if indel.is_deletion():
        start = indel.position
        end = indel.position + len(indel.reference) - 1
        seq = seq[:start] + seq[end:]
    else:
        seq = seq[:indel.position] + \
            list(indel.alteration[1:]) + seq[indel.position:]
    return seq


def select_snps(vcf: dict[int, list[Mutation]], lims: tuple[float, float]) -> list[Mutation]:
    snps = []
    for pos in vcf:
        snp = select_snp(vcf[pos], lims)
        snps.append(snp)
    return snps


def select_snp(muts: list[Mutation], lims: tuple[float, float]) -> Mutation:
    min_freq, max_freq = lims
    snps = [m for m in muts if not m.is_indel()]
    snps.sort(key=lambda m: m.frequency, reverse=True)

    if not snps:
        return None

    if snps[0].frequency >= max_freq:
        return snps[0].alteration

    ref_freq = 100
    degs = ''

    for snp in snps:
        if snp.frequency >= min_freq:
            degs += snp.alteration
        ref_freq -= snp.frequency
    if ref_freq >= max_freq:
        return snp.reference
    if ref_freq >= min_freq:
        degs += snp.reference
    return to_iupac(degs)


def degenerate_snp(snps: list[Mutation], max_freq: float) -> Mutation:
    deg_str = ''
    for snp in snps:
        if snp.frequency >= max_freq:
            return snp


def to_iupac(nucls: str) -> str:
    sorted_nucls = ''.join(sorted(nucls))
    return iupac_names[sorted_nucls]


iupac_names = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y',  'GT': 'K',
    'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B',
    'ACGT': 'N'}


if __name__ == '__main__':
    cli()
