---
aliases:
  - Methylation microarrays
url: https://www.illumina.com/techniques/microarrays/methylation-arrays.html
---

Los **arrays de metilación** proporcionan información de la metilación en sitios CpG específicos (también algunos no CpG) de todo el genoma (genome-wide). Los arrays más populares son los fabricados por Illumina, que se diferencian en el número de sitios CpG interrogados:

- Infinium HumanMethylation450K BeadChip (popular en datos del TCGA)
- Infinium HumanMethylationEPIC BeadChip (850K)
- [Infinium HumanMethylationEPIC v2.0 BeadChip](https://emea.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html) (900K)

La química de los arrays de Infinium depende de una [[Conversión por bisulfito]] previa de la muestra. Tras una amplificación del genoma completo, las citosinas no metiladas son convertidas a timinas con un intermediario de uracilo.

Las sondas presentan regiones de 50 bases complementarias a la secuencia convertida por bisulfito, con el sitio CpG que se va ensayar en el extremo 3' de la secuencia. Tras la hibridación con la secuencia de DNA convertida por bisulfito, se hace una extensión de base única en presencia de terminadores ddNTP marcados con fluorescencia. Los principales arrays mencionados combinan dos tipos de sondas:

- Infinium I. Por cada sitio CpG hay dos secuencias de sonda (una para la secuencia metilada y otra para la secuencia sin metilar). Ocupan el doble de espacio físico pero son útiles en regiones densas en CpG.
- Infinium II. Por cada sitio CpG hay una sola secuencia de sonda.

![[Pasted image 20240325191901.png]]

En todo caso, las sondas y los terminadores están diseñados para que la extensión de una base en locus metilados de lugar a señal verde, mientras que la extensión de una base en locus no metilados da lugar a señal roja. De esta forma, la relación rojo/verde en una posición será indicativa del grado de metilación del sitio en una muestra. 

Los datos crudos de intensidad de señal vienen en formato [[IDAT format|IDAT]].

## References

- [Infinium Human Methylation BeadChip](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/ewas-suite/tutorial.html)
- [Illumina HumanMethylation BeadChip for Genome-Wide DNA Methylation Profiling: Advantages and Limitations | SpringerLink](https://link.springer.com/referenceworkentry/10.1007/978-3-319-31143-2_89-1)
- [Validation of the new EPIC DNA methylation microarray (900K EPIC v2) for high-throughput profiling of the human DNA methylome - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9988339/)
- [Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/)