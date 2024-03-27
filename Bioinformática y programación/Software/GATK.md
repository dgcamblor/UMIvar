---
url:
---
GATK has good documentation on best practices.

[Best Practices Workflows – GATK](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)

## Mutect2

The variant calling process that Mutect2 uses employ assembled haplotypes ( #write ). To visualize the BAM specifically created by Mutect2, you can use the option `--bam-output`.

- `--dont-use-soft-clipped-bases`. It is useful in amplicon sequencing, for example, were sequencing is restricted to very specific regions.

## FilterMutectCalls

Allows for filtering the raw variant calling output of mutect2.

- `--max-events-in-region`. The default of `2` may filter too harshly. It can be adjusted to more.

## Picard

Starting with version 4.0, GATK contains a copy of the [Picard](http://broadinstitute.github.io/picard/) toolkit, so all Picard tools are available from within GATK itself.

![[Picard]]

