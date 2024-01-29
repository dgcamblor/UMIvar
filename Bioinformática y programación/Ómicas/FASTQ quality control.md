
FASTQ files will be encountered in any NGS application (genomics, transcriptomics, epigenomics...). Quality of base calling is encoded in the Phred score (line 4 of FASTQ files):

$$
Q = -10 \log_{10} P
$$

where $P$ is the probability of the base being called incorrectly. A guide to Phred scores is:

| Phred score | Probability of incorrect base call | Base call accuracy |
| ----------- | --------------------------------- | ------------------ |
| 10          | 1 in 10                            | 90%                |
| 20          | 1 in 100                           | 99%                |
| 30          | 1 in 1000                          | 99.9%              |
| 40          | 1 in 10000                         | 99.99%             |

## Software for FASTQ QC

- [[FastQC]]
- Fastq-screen
- [[MultiQC]]

## References

- [Quality control using FastQC | Training-modules](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html)