Unique Molecular Identifiers (UMIs) are short nucleotide sequences added to a sample prior to the amplification (preparation of the sequencing library). 

## Bioinformatic processing

There are two main bioinformatic approaches for using UMI (Bohers et al., 2021):

- **Grouping** and **deduplicating** the duplicate reads (those that have the same UMI). Tools for that purpose are [[UMI-tools]] and [[fgbio]].
- **Variant calling** that takes into account the information of the UMI duplicates, without deduplication.

## References
Bohers, Elodie, Pierre-Julien Viailly, and Fabrice Jardin. “cfDNA Sequencing: Technological Approaches and Bioinformatic Issues.” _Pharmaceuticals_ 14, no. 6 (June 2021): 596. [https://doi.org/10.3390/ph14060596](https://doi.org/10.3390/ph14060596).