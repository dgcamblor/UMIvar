---
aliases:
  - Methylation microarrays
url: https://www.illumina.com/techniques/microarrays/methylation-arrays.html
---

Los **arrays de metilación** proporcionan información de la metilación en sitios CpG específicos (también algunos no CpG) de todo el genoma (genome-wide). Los arrays más populares son los fabricados por Illumina, que se diferencian en el número de sitios CpG interrogados:

- Infinium HumanMethylation450K
- [Infinium MethylationEPIC v2.0 Kit](https://emea.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html)

La química de los arrays de Infinium depende de una [[Conversión por bisulfito]] previa de la muestra. Tras una amplificación del genoma completo, las citosinas no metiladas son convertidas a timinas con un intermediario de uracilo.

En el chip, cada posición es interrogada por dos tipos de sondas que tienen ambas una secuencia de 50 nucleótidos con: 1) un tramo común complementario a la secuencia previa al sitio CpG; 2) un nucleótido final complementario o bien a la secuencia metilada tras secuenciación de bisulfito (una T) o bien a la secuencia no metilada (una C). De esta forma:

- Methylation (M) leads to red fluorescent signal.
- Unmethylation (U) leads to green fluorescent signal.

![[Pasted image 20231129144129.png]]

The output is an [[IDAT format]] file.

## References

- [Infinium Human Methylation BeadChip](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/ewas-suite/tutorial.html)