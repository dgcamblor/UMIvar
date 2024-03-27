---
aliases:
  - Methylation microarrays
url: https://www.illumina.com/techniques/microarrays/methylation-arrays.html
---

Los **arrays de metilación** proporcionan información de la metilación en sitios CpG específicos (también algunos no CpG) de todo el genoma (genome-wide). Los arrays más populares son los fabricados por Illumina, que se diferencian en el número de sitios CpG interrogados:

- Infinium HumanMethylation450K
- [Infinium MethylationEPIC v2.0 Kit](https://emea.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html)

La química de los arrays de Infinium depende de una [[Conversión por bisulfito]] previa de la muestra. Tras una amplificación del genoma completo, las citosinas no metiladas son convertidas a timinas con un intermediario de uracilo.

En el chip, cada posición es interrogada por dos tipos de sondas que tienen ambas una secuencia de 50 nucleótidos con: 1) un tramo común complementario a la secuencia previa al sitio CpG; 2) un nucleótido final complementario o bien a la secuencia metilada tras secuenciación de bisulfito (una T) o bien a la secuencia no metilada (una C). La muestra de ADN tratada con bisulfito se hace hibridar con la sondas en los pocillos del chip. Se realiza una extensión de cadena única en la que las sondas actúan como cebadores para una reacción de extensión de una base utilizando terminadores de cadena marcados con fluorescente. Los terminadores T y A tienen fluorescentes rojos, y los terminadores G y C tienen fluorescentes verdes. De esta forma, en cada sitio:

- La metilación (M) lleva a señal fluorescente roja.
- La no methylation (U, unmethylation) lleva a señal fluorescente verde.

![[Pasted image 20231129144129.png]]

La relación rojo/verde será indicativa del grado de metilación del sitio en una muestra. Los datos crudos de intensidad de señal vienen en [[IDAT format]].

The output is an [[IDAT format]] file.

## References

- [Infinium Human Methylation BeadChip](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/ewas-suite/tutorial.html)