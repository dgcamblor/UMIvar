# UMIvar: UMI-aware simulation of variants in BAM files

![Version](https://img.shields.io/github/v/tag/dgcamblor/UMIvar?label=Version)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<img src="https://github.com/dgcamblor/UMIvar/assets/74027422/53e0c724-3678-4710-9056-060dc37b056b" width=50% height=50%>

## Description

**UMIvar** is a tool for simulating mutations (SNVs, MNVs and indels) in BAM files in a UMI-aware manner. UMIs are used in NGS to identify unique reads among PCR duplicates. When a specific read reports a variant to the reference genome, it is expected that all the PCR duplicates of that read will also report the same variant. This is not always the case when the variant is introduced artificially, as it can occur in different steps of the amplification process or due to punctual errors in the sequencing. In these cases, the variant may not be present in all the PCR duplicates of the read, and end up in a subset of the PCR duplicates. UMIvar takes this into account and simulates the variant in all the PCR duplicates of the read.

UMIvar is specifically designed to be used for the benchmarking of variant calling pipelines.

## Installation

UMIvar is written in Python 3 (works for Python 3.10+) and requires the following python packages:

- pysam
- numpy
- umi_tools (only for the `--edit_threshold` parameter)
- pyfaidx (only if using random variants)

Installation of these dependencies can be easily done using pip:

```bash
pip3 install -r requirements.txt
```

It also requires the [samtools](http://www.htslib.org/) executable to be in the PATH. Check the installation instructions [here](http://www.htslib.org/download/).

## Input

UMIvar takes as an input an NGS sample in the form of a **BAM file**, which must have its index (`.bai` file) in the same folder. This file must contain UMI tagged reads, with the UMI tag included at the end of each query name, separated from the rest of the metainformation by a non letter character (e.g. "_" or ":"). For example, a read with the following query name would be valid:

```text
@MISEQ:177:C3Y1KACXX:2:1101:1465:1937_GATTAGATT
```

On the other hand, the program reads the variants to be introduced from a **CSV** or a **VCF file**. If variants are stored in a CSV file, it must contain the following columns:

```text
chromosome,position,reference,alternative,allele_frequency
```

If variants are stored in a VCF file, it must contain the usual fields, and the allele frequency must be included in the `AF` tag under the `INFO` field.

Note that UMIvar supports:

- SNVs, MNVs and indels.
- Several variants in the same read.

It is recommended that variants are sorted by chromosome and position, and they should be expressed in accordance with [normalization rules](https://genome.sph.umich.edu/wiki/Variant_Normalization).

Alternatively, random variants can be introduced in the BAM file by specifying the number of variants to be introduced. In this case, a reference genome must be provided in FASTA format, and the variants will be randomly generated.

## Usage

UMIvar can be launched from the command line using the following command:

```bash
python3 bin/umivar.py -i <INPUT_BAM> -v <VARIANTS> [OPTIONS]
```

The required arguments are:

- `-i` or `--input_bam`: the input BAM file (with the aforementioned specifications) where the variants will be introduced.
- `-v` or `--variants`: the CSV or VCF file with the variants to be introduced, following the aforementioned specifications. Alternatively, the number of random variants to be introduced can be specified here: inputting a number will activate the random variant mode (requires `-f`).

The optional arguments are:

- `-o` or `--output_bam`: the path to the output file. **Default:** `<input_bam_name>_UV.bam`.
- `-t` or `--output_type`: the format of the output, to be chosen from `BAM`, `SORTED_BAM` or `FASTQ`. **Default:** `BAM`.
- `-b` or `--bed`: a BED file with the regions to be considered for the simulation. The output BAM file will only contain reads that map to these regions. **NOTE:** A hash of the input files is generated and saved to avoid reprocessing the same files. **Default:** all the reads in the input BAM file are considered.
- `-f` or `--ref_fasta`: the path to the reference genome in FASTA format. Required only if random variants are to be introduced.
- `-e` or `--edit_threshold`: Edit distance threshold for the UMI clustering using UMI-tools (accounting for errors in the UMI sequence). If `0`, no clustering is performed. **Default:** `0`.
- `-s` or `--seed`: the seed for the random number generator. **Default:** `12`.

## Output

UMIvar will examine the reads covering each variant and will introduce those variants at the nearest possible allele frequency to the one specified in the CSV file, depending on the available UMIs (or UMI groups). The UMI simulation tries to reduce the strand bias of the introduced variants. Finally, it will write the resulting reads in the specified output BAM file (or FASTQ).

As the final allele frequency may differ from the one specified in the CSV file, and as some of the variants may not be covered enough by the reads and discarded, a CSV file with the information of the variants that were finally introduced is also generated and saved to the output folder (`<output_bam_name>.csv`). This file contains following columns:

```text
chromosome,position,reference,alternative,achieved_af,var_coverage,umis,strand_bias
```

Recording the achieved allele frequency (`achieved_af`), the number of reads covering the variant (`var_coverage`), the number of UMIs in which the variant was included (`umis`) and the strand bias (`strand_bias`). This information can help to assess the performance of the variant calling pipeline.

UMIvar can generate two FASTQ files with the paired reads resulting from the simulation, ready to be realigned to the reference genome (`<output_bam_name>_R1.fastq` and `<output_bam_name>_R2.fastq`).

## Infographic

![UMIvar infographic](https://github.com/dgcamblor/UMIvar/assets/74027422/b3c6c3b7-3631-4cb0-a298-71a0105d0c87)

## License

Distributed under the MIT License (for more information, see the `LICENSE` file).