---
githubUrl: https://github.com/Daniel-Liu-c0deb0t/UMICollapse
paperUrl: https://peerj.com/articles/8275/
dedupStrategy: Consensus
---
UMI
## Installation

```
# Dependencies - Java 11
mkdir lib
cd lib
curl -O -L https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar
curl -O -L https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
cd ..

# UMICollapse
conda install -c bioconda umicollapse
```

## Input

A BAM file.
## Usage

```
./umicollapse bam -i paired_example.bam -o dedup_paired_example.bam --umi-sep : --paired --two-pass
```

Which is the equivalent of [[UMI-tools]] `dedup`.

> [!Note]
> Fast and memory efficient.