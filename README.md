# UMIvar: UMI-aware introduction of simulated variants in BAM files

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description

**UMIvar** is a tool for simulating mutations (SNVs and indels) in BAM files in a UMI-aware manner. UMIs are used in NGS to identify unique reads among PCR duplicates. When a specific read reports a variant to the reference genome, it is expected that all the PCR duplicates of that read will also report the same variant. This is not always the case when the variant is introduced artificially, as it can occur in different steps of the amplification process and end up in a subset of the PCR duplicates. UMIvar takes this into account and simulates the variant in all the PCR duplicates of the read.

UMIvar is specifically designed to be used for the benchmarking of variant calling pipelines.

## Installation

UMIvar is written in Python 3 (works for Python 3.10+) and requires the following python packages:

- pysam
- numpy

Installation of these dependencies can be easily done using pip:

```bash
pip3 install -r requirements.txt
```

It also requires the [samtools](http://www.htslib.org/) executable to be in the PATH. Check the installation instructions [here](http://www.htslib.org/download/).

## Input

UMIvar takes as an input an NGS sample in the form of a **BAM file**, which must have its index (`.bai` file) in the same folder. This file must contain UMI tagged reads, with the UMI tag included at the end of each query name, separated from the rest of the metainformation by a non letter character (e.g. "_" or ":"). For example, a read with the following query name:

```text
@MISEQ:177:C3Y1KACXX:2:1101:1465:1937_GATTAGATT
```

On the other hand, the program reads the variants to be introduced from a **CSV file**. The CSV file must contain the following columns:

```text
chromosome,position,reference,alternative,allele_frequency
```

Variants should be sorted by chromosome and position, and they should be expressed in accordance with [normalization rules](https://genome.sph.umich.edu/wiki/Variant_Normalization).

## Usage

UMIvar can be launched from the command line using the following command:

```bash
python3 bin/umivar.py -i <INPUT_BAM> -v <VARIANTS> [OPTIONS]
```

The required arguments are:

- `-i` or `--input_bam`: the input BAM file with the aforementioned specifications.
- `-v` or `--variants`: the CSV file with the variants to be introduced, following the aforementioned specifications.

The optional arguments are:

- `-o` or `--output_bam`: the output BAM file. **Default:** `<input_bam_name>_UV.bam`
- `-b` or `--bed`: a BED file with the regions to be considered for the simulation. The output BAM file will only contain reads that map to these regions. **Default:** all the reads in the input BAM file are considered
- `-f` or `--freq_mode`: the mode for the allele frequency (both input and output). **Default: dedup**. The possible values are:
  - `dedup`: the allele frequency is computed after removing PCR duplicates (reads with equal UMIs).
  - `dup`: the allele frequency is considered before removing PCR duplicates. Allows to introduce variants in the BAM as if it had no UMI tags, while still keeping the UMI information for the simulation.
- `-s` or `--seed`: the seed for the random number generator. **Default:** `12`.
- `-c` or `sb_correction`: perform a simple strand bias correction whenever possible. **Default:** `False`.

## Output

UMIvar will examine the reads covering each variant and will introduce those variants at the nearest possible allele frequency to the one specified in the CSV file. It will then write the resulting reads in the specified output BAM file.

As the final allele frequency may differ from the one specified in the CSV file (it is subject to the available UMIs and their frequencies), and as some of the variants may not be covered enough by the reads and discarded, a CSV file with the information of the variants introduced is also generated and saved to the output folder (`<output_bam_name>.csv`). This file contains the same columns as the input CSV file, in addition to the number of umis in which the variant was included and the strand bias:

```text
chromosome,position,reference,alternative,achieved_af,umis,strand_bias
```

UMIvar will also generate two FASTQ files with the paired reads resulting from the simulation, ready to be realigned to the reference genome (`<output_bam_name>_R1.fastq` and `<output_bam_name>_R2.fastq`).

## License

Distributed under the MIT License (for more information, see the `LICENSE` file).
