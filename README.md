# Consenser

Create a consensus sequence from a reference FASTA and a LoFreq-style VCF.

## Installation

```bash
git clone https://github.com/izsvenezie-virology/consenser.git
cd consenser
python3 -m pip install .
```

## Usage

```bash
consenser [OPTIONS] REFERENCE_FILE VCF_FILE
```

## Options

| Option                             | Description                                                                            | Default  |
| ---------------------------------- | -------------------------------------------------------------------------------------- | -------- |
| `-o`, `--output FILE`              | Output FASTA file. Use `-` for stdout.                                                 | `-`      |
| `-c`, `--coverage-file FILE`       | Tab-separated coverage file with columns: `Chrom`, `Position`, `Coverage` (no header). | not used |
| `-H`, `--header TEXT`              | Output sequence header template. `CHROM` is replaced with the original contig name.    | `CHROM`  |
| `-m`, `--minimum-coverage INTEGER` | Minimum coverage required to keep a base; lower values are masked as `N`.              | `10`     |
| `-s`, `--snp-threshold DECIMAL`    | Minimum SNP allele frequency to include in IUPAC consensus.                            | `0.25`   |
| `-i`, `--indel-threshold DECIMAL`  | Minimum INDEL allele frequency to apply insertion/deletion.                            | `0.50`   |
| `-v`, `--version`                  | Show version information.                                                              | -        |
| `-h`, `--help`                     | Show help.                                                                             | -        |

## Inputs

- `REFERENCE_FILE`: FASTA reference used during alignment.
- `VCF_FILE`: VCF file (expected LoFreq-style `AF=` entries in `INFO`).

## Example

```bash
consenser \
  --coverage-file sample.cov \
  --minimum-coverage 20 \
  --snp-threshold 0.30 \
  --indel-threshold 0.60 \
  --header "sample_CHROM" \
  --output consensus.fasta \
  reference.fasta variants.vcf
```

## Degenerations

At each position, variants are filtered by thresholds:

- SNPs with `AF < --snp-threshold` are excluded from the IUPAC set.
- INDELs with `AF < --indel-threshold` are ignored.

For SNP calling, `consenser` then builds the consensus base as follows:

- If one or more SNPs (including reference) pass the SNP threshold, they are converted to an IUPAC code.
- If no SNP passes the threshold, the most frequent SNP is used.
- If two or more SNPs are tied as most frequent, a degenerate IUPAC base is still produced.

### Avoid Degenerations

Set `-s`, `--snp-threshold` to `1` to force selection of the most frequent SNP at each position.
Degenerations can still appear when two or more SNPs have exactly the same highest frequency.
