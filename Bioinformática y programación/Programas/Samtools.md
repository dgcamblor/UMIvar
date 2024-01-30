
## [view](https://www.htslib.org/doc/samtools-view.html)

`samtools view` is used to read a file in format SAM/BAM.

```
samtools view <file>
```

- `-b`: Output in BAM format. Useful for conversion SAM -> BAM (`samtools view -b <SAM> > <BAM>)
- `-h`: Output with BAM header ([[BAM format#Header]]).
- `-q`: Do not output reads with [[Mapping quality (MAPQ)]] less than the specified.

### SAM to BAM

```bash 
samtools view -bS SAMPLE.sam > SAMPLE.bam
```

### BAM to SAM

```bash
samtools view -h SAMPLE.bam > SAMPLE.sam
```

## fastq

### BAM to fastq

First, the BAM file must be sorted by name.

```bash
samtools sort -n SAMPLE.bam -o SAMPLE.nsorted.bam
```

Then, the FASTQ files can be generated.

```bash
samtools fastq -@ 8 SAMPLE.nsorted.bam \
-1 SAMPLE_R1.fastq.gz \
-2 SAMPLE_R2.fastq.gz \
-0 /dev/null -s /dev/null -n
```

