
## concat

This command is used to concatenate or combine VCF/BCF files that **belong to the same sample or samples**, they should have the same sample columns appearing in the same order. It is suitable for combining files split by region or chromosome, i.e., the same samples.

## merge

This command is used to put together files that come from **different samples or different sets of samples**.

## norm (variant normalization)

```bash
bcftools norm -f ref.fa -O z -m any variants.vcf
```

- `-m`: Normalization method. `any` allows any indel to be left-aligned (even if it is not parsimonious).
- `-f`: Path to FASTA file.
- `-O`: Output type (`z` = gzip).

