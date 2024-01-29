
## [view](https://www.htslib.org/doc/samtools-view.html)

`samtools view` is used to read a file in format SAM/BAM.

```
samtools view <file>
```

- `-b`: Output in BAM format. Useful for conversion SAM -> BAM (`samtools view -b <SAM> > <BAM>)
- `-h`: Output with BAM header ([[BAM format#Header]].
- `-q`: Do not output reads with [[Mapping quality (MAPQ)]] less than the specified.


## fastq

Converting BAM to FASTQ

```bash
samtools sort -n SAMPLE.bam -o SAMPLE.nsorted.bam
```

```bash
samtools fastq -@ 8 SAMPLE.nsorted.bam \
-1 SAMPLE_R1.fastq.gz \
-2 SAMPLE_R2.fastq.gz \
-0 /dev/null -s /dev/null -n
```

