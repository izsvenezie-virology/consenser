import click
from click.types import File

@click.command()
@click.option('-c', '--coverage', type=File('r'), required=True)
@click.option('-r', '--reference', type=File('r'), required=True)
@click.option('-v', '--vcf', type=File('r'), required=True)
def cli(coverage, reference, vcf):
    cov_by_chrom = read_coverage(coverage)
    ref_by_chrom = read_reference(reference)
    vcf_by_chrom = read_vcf(vcf)
    for chrom in vcf_by_chrom:
        for pos in vcf_by_chrom[chrom]:
            for mut in vcf_by_chrom[chrom][pos]:
                print(f'Pos: {pos}; Mut: {mut}')


class Mutation():
    def __init__(self, pos, ref, alt, freq) -> None:
        self.position = int(pos)
        self.reference = ref
        self.alteration = alt
        self.frequency = float(freq)*100
    
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
        infos = info.split(';')
        af = [i for i in infos if i.startswith('AF=')][0][3:]
        mut = Mutation(pos, ref,alt, af)

        if chrom not in chroms_vcf:
            chroms_vcf[chrom] = {}
        if pos not in chroms_vcf[chrom]:
            chroms_vcf[chrom][pos] = [mut]
        else:
            chroms_vcf[chrom][pos].append(mut)
    return chroms_vcf


if __name__ == '__main__':
    cli()
