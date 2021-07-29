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
def cli(coverage, reference, vcf, deg, no_deg, min_cov, width, split, alter_names):
    if not deg and not no_deg:
        raise ClickException("'-d'/'--deg' or '-n'/--no-deg' is required")
    if deg and no_deg:
        raise ClickException("Use '-d'/'--deg' or '-n'/--no-deg', not both")
    if (no_deg):
        deg = (no_deg, no_deg)

    cov_by_chrom = read_coverage(coverage)
    ref_by_chrom = read_reference(reference)
    vcf_by_chrom = read_vcf(vcf)
    chroms = ref_by_chrom.keys()

    for chrom in chroms:
        cons_seq = create_consensus_sequence(
            ref_by_chrom[chrom], vcf_by_chrom[chrom], cov_by_chrom[chrom], deg, min_cov)


class Mutation():
    def __init__(self, chrom, pos, ref, alt, freq) -> None:
        self.chromosome = chrom
        self.position = pos
        self.reference = ref
        self.alteration = alt
        self.frequency = freq

    def is_indel(self):
        return len(self.reference) > 1 or len(self.alteration) > 1

    def is_deletion(self):
        return len(self.reference) > 1

    def is_insertion(self):
        return len(self.alteration) > 1

    def __str__(self) -> str:
        return(f'Pos: {self.position}; Ref: {self.reference}; Alt: {self.alteration}; Freq: {self.frequency}')


def read_coverage(cov_file: File) -> dict[str, dict[int, int]]:
    chroms_cov = {}
    for line in cov_file.read().split('\n'):
        if line == '':
            continue
        chrom, pos, cov = line.split('\t')
        if chrom not in chroms_cov:
            chroms_cov[chrom] = {}
        chroms_cov[chrom][int(pos)] = int(cov)
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


def read_vcf(vcf_file: File) -> dict[str, dict[int, list[Mutation]]]:
    chroms_vcf = {}
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


def create_consensus_sequence(seq: str, vcf: dict[int, list[Mutation]], cov: dict[int, int], lims: tuple[float, float], min_cov: int):
    min_freq = lims[0]
    max_freq = lims[1]
    min_indel = 50.0

    seq = list(seq)

    for pos in range(len(seq), 0, -1):
        index = pos - 1
        if cov[pos] < min_cov:
            seq[index] = 'N'
            continue
        if pos in vcf:
            muts = vcf[pos]
            check_reference(seq, muts)
            indel = get_consensus_indel(muts, min_indel)
            if indel:
                seq = apply_indel(seq, indel)
            snps = get_consensus_snps(muts, min_freq)
            if snps:
                snp = degenerate_snp(snps, max_freq)
    return ''.join(seq)


def check_reference(seq: list[str], muts: list[Mutation]) -> None:
    '''Checks if reference and mutation's reference match'''
    index = muts[0].position - 1
    if seq[index] != muts[0].reference[0]:
        msg = f'Mismatch between mutation and reference at position {muts[0].position}'
        raise ValueError(msg)


def get_consensus_indel(muts: list[Mutation], freq_min: float) -> Mutation:
    '''Returns the indel to insert/delete from consensus sequence'''
    indels = [m for m in muts if m.is_indel() and m.frequency >= freq_min]
    indels.sort(key=lambda m: m.frequency, reverse=True)
    return indels[0] if indels else None


def get_consensus_snps(muts: list[Mutation], freq_min: float) -> list[Mutation]:
    '''Returns the mutations over the minimum frequency'''
    return [m for m in muts if not m.is_indel() and m.frequency >= freq_min]


def apply_indel(seq: list[str], indel: Mutation) -> list[str]:
    if indel.is_deletion():
        start = indel.position
        end = indel.position + len(indel.reference) - 1
        seq = seq[:start] + seq[end:]
    else:
        seq = seq[:indel.position] + list(indel.alteration[1:]) + seq[indel.position:]
    return seq


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
