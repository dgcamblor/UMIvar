---
tags:
  - program
---

```bash
bwa mem
```

- `-t`: Number of threads to use for the alignment.
- `-p`: Smart pairing. When two reads are adjacent in the same FASTQ file and they have the same same, they will be considered a read pair. Useful in pipe commands. Not much if paired reads are properly separated.

- `-K`: Reads are processed in chunks, regardless of the number of threads, making `bwa` deterministic and reproducible with different thread counts.
- `-M`: Marks shorter split hits as secondary alignments instead of supplementary alignments ([[Chimeric alignment]]). Useful only for Picard compatibility (not needed in recent versions).
- `-Y`: Use [[Soft-clipped]] bases in supplemental alignments ([[Chimeric alignment]]) instead of [[Hard-clipped]] bases (the default output).
- `-R`: Sets the **read group** header line. Tab separated fields. This labels the read groups for better downstream processing. See: [[Read group]].