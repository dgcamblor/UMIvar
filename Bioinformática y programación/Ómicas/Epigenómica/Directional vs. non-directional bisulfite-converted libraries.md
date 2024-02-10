
Two main types of libraries can be generated in bisulfite sequencing:

- Directional. The most followed approach.
- Non-directional

## Directional bisulfite sequencing

In directional bisulfite sequencing, only the original top and bottom strands of DNA are sequenced. It is the most common approach, and the one most aligners ([[Bismark]], [[bwa-meth]]) are adapted for.

The forward sequencing reads will correspond to a bisulfite-converted version of either the original top or the original bottom strand (the C-to-T reads) and the reverse sequencing reads will correspond to the complement of the original top or the complement of the original bottom strand (the G-to-A reads).

### Example: top strand

For example, if the original top strand is:

```{r}
5' - ATCGATCG - 3'
```

After bisulfite conversion, the unmethylated C's will be converted to T's, and the original top strand (forward read) will be:

```{r}
5' - ATTGATTG - 3'
```

With its complement (reverse read) being:

```{r}
3' - TAACTAAC - 5'
```

### Example: bottom strand

Following the same example, if the original bottom strand is:

```{r}
3' - TAGCTAGC - 5'
```

After bisulfite conversion, the unmethylated C's will be converted to T's, and the original bottom strand (forward read) will be:

```{r}
3' - TAGTTAGT - 5'
```

With its complement (reverse read) being:

```{r}
5' - ATCAATCA - 3'
```

## Non-directional bisulfite sequencing

Non-directional bisulfite sequencing, on the other hand, does not preserve the information about which strand a read comes from. The methylation status of cytosines on both strands is determined together. 