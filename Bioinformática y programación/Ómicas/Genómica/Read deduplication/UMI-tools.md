---
githubUrl: https://github.com/CGATOxford/UMI-tools
paperUrl: http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract
dedupStrategy: Highest MAPQ
tags:
  - program
---

## Installation

```
conda install -c bioconda umi_tools
```


A quick guide on using UMI-tools can be found at: https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md

## Input

A BAM file.

## Usage
### extract

Extracting the UMI from the sequence read and adding it to the read name.

```bash
umi_tools extract --stdin=reads.fastq.gz --bc-pattern=NNNNNNNN --stdout=processed.fastq.gz --log=processed.log 
```

### dedup

Deduplicating the reads based on the UMI.

```bash
umi_tools dedup -I example.bam -S deduplicated.bam --output-stats=deduplicated
```

### group (optional)

Grouping the reads based on the UMI can optionally be done before deduplication. It essentially has two usages:

- To output a tagged BAM file with the groups
- To output a flatfile that describes read groups (i.e. how many duplicates).

```bash
umi_tools group -I mapped.bam --paired --output-bam -S mapped_grouped.bam --group-out=groups.tsv
```

If `--output-bam` is not specified, only the stats will be outputted. This can be useful if you want to keep track of the number of reads per UMI.