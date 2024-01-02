---
githubUrl: https://github.com/fulcrumgenomics/fgbio
paperUrl: Not published
dedupStrategy: Consensus
tags:
  - program
---

## Installation

```bash
wget https://github.com/fulcrumgenomics/fgbio/releases/download/2.1.0/fgbio-2.1.0.jar
mkdir fgbio
mv fgbio-2.1.0.jar fgbio/
```

## Input

For the best practice pipeline, the input should be a `FASTQ` with UMIs, specifying the read structure with `--read-structures`. UMIs can also be extracted from read names using `--extract-umis-from-read-names` in the `FastqToBam` step.

UMI appended at the end, separated by `:`, 8th field.
## Best practice consensus pipeline

[best-practice-consensus-pipeline.md](https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md)

- **FASTQ -> Grouped BAM.** 
- **Grouped BAM -> Filtered Consensus**