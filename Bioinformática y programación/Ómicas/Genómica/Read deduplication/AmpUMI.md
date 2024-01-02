---
githubUrl: https://github.com/pinellolab/AmpUMI
paperUrl: https://academic.oup.com/bioinformatics/article/34/13/i202/5045706?login=false
dedupStrategy: Most frequent sequence
---

## Installation

```bash
pip install sympy mpmath numpy  # Dependencies
pip install git+https://github.com/pinellolab/AmpUMI.git
```

## Input

Single FASTQ with R1 and R2 merged.

> First, merge your R1 and R2 using [FLASh](https://ccb.jhu.edu/software/FLASH/) or a similar tool. These merged reads can then be used as input for AmpUMI.
> For scRNA-seq or other applications where paired reads may not overlap (and cannot be merged before deduplicating) there are other tools (e.g. [UMI-tools](https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md#paired-end-sequencing)) that are designed for those use cases.

This renders the software not useful for the [[Comparative analysis of UMI deduplication software in variant calling error correction]].
## Use

Using the process mode with `python AmpUMI.py`. The Process mode will parse and trim the UMI. Performs two error correction steps:

1. The most frequent amplicon sequence for each UMI will be accepted as the consensus sequence and the other sequences will be assumed to come from sequencing error and will be discarded. 
2. (Optional) Correction of errors in the UMI.

```bash
python3 AmpUMI --fastq input.fastq --umi_regex "^IIIII" --fastq_out input.fastq.dedup.fastq
```

- `--fastq` is the input fastq file
- `--umi_regex` is the regular expression to match the UMI (in this case, the first 5 bases of the read)
- `--fastq_out` is the output fastq file