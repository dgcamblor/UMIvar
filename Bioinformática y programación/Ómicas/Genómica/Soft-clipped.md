**Soft-clipped** bases are bases (in the 5' and 3' extremes of a read) that are not part of an alignment, although they are included in the read sequence (contrary to [[Hard-clipped]] bases).

A soft-clipped base can mean something of biological origin for which the aligner simply couldn't assign a location. Given that definition, it is actually useful to not remove said bases, and take them into account for variant calling processes in some cases (they may be due to indels).

## References

- [Dave's Wiki | SAM](https://davetang.org/wiki/tiki-index.php?page=SAM)