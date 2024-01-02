
## norm (variant normalization)

```bash
bcftools norm -f ref.fa -O z -m any variants.vcf
```

- `-m`: Normalization method. `any` allows any indel to be left-aligned (even if it is not parsimonious).
- `-f`: Path to FASTA file.
- `-O`: Output type (`z` = gzip).

