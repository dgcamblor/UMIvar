---
tags:
  - phd
---

## Functions

`fastq_preprocessing()`

- {1} R1
- {2} R2
- {3} sample
- {4} tmp_dir
- {5} analysis_dir
- {6} threads
- {7} has_umi
- {8} umi

`mapping_and_bampostprocessing()`

- {1} R1.fastq = `${tmp_dir}/QC/filtered_${sample}_1.fastq.gz`
- {2} R2.fastq =  `${tmp_dir}/QC/filtered_${sample}_2.fastq.gz`
- {3} tmp_dir
- {4} analysis_dir
- {5} threads
- {6} sample
- {7} genome
- {8} panel (panel folder)
- {9} has_umi
- {10} amplicon
- **NEW:** {11} umi

## fgbio implementation

Take care of `-M` flag (`mark shorter split hits as secondary`). This can cause problems with [[fgbio]] tools. Mapping command should be:

```bash
bwa mem -t 16 -p -K 150000000 -Y ref.fa
```

Approximations:
 
 - **FastqToBam** (without --read-structures) -> **AnnotateBamWithUmis** (needs additional {10} variable = umi)
 -  Extracting UMIs with UMI-tools (normal pipeline) and replacing `_` for `:`. Script needed. Using **FastqToBam** (--extract-umis-from-read-names).
 - Extracting UMIs with **bcl2fastq** (uses ":") -> Using **FastqToBam** (--extract-umis-from-read-names).
 - Introducing the UMI in the sequence of the reads with: `join <(zcat EN104_UMI.fastq.gz | nl ) <(zcat EN104_1.fastq.gz | nl) | awk -F ' ' '{if ($2 ~ /^@ST-E00114/ ) {print $4" "$5;} else {print $2$3;}}' | gzip > EN104_1.umi.fastq.gz`