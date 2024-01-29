
Read deduplication is the process of correcting for [[Read duplicates]].

Whenever UMIs are available, UMI deduplication software is needed:

- [[UMI-tools]] -> Transcriptomics
- [[fgbio]] -> Variant calling

When **not working with amplicon sequencing**, MarkDuplicates can be used.

In **amplicon sequencing**, in those cases where UMIs are not available, deduplication is not performed because it's impossible to distinguish between genuine reads from different cells and PCR duplicates. Deduplication based solely on mapping position and orientation would only allow one alignment to a given position.