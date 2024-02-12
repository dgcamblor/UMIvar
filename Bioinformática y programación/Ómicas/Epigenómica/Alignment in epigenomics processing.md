Alignment to the reference genome depends on the strategy used to determine methylation.

## Bisulfite conversion-based methods

Requires special aligners, as the sequence itself has been altered by the bisulfite conversion and does not align directly to the reference genome. Two algorithms are available: wild card and three-letter.

- Wild card algorithm. Both Cs and Ts map into Cs in the reference genome.
- Three-letter algorithm. Converts all Cs in the reference genome and the reads into Ts, so that standard aligners can be used.

Three-letter aligners seem to outperform wildcard aligners in running time and peak memory usage [@gongAnalysisPerformanceAssessment2022]. The most popular letter algorithms are:

- [[Bismark]]. Extensively used software for the analysis of bisulfite sequencing data, which includes alignment, duplicate removal (requires [[Picard]], samtools), and methylation calling.
- BS-Seeker2
- BWA-Meth

> [!NOTE] EpiQC study insights
> The SEQC2 study (and additional studies) seem to suggest that bwa-meth has the highest uniquely mapped read rates and the lowest unmapped reads [@fooxSEQC2EpigenomicsQuality2021; @gongAnalysisPerformanceAssessment2022].

## Enrichment-based methods

The standard aligners can be used.