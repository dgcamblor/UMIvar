## INFO

Contiene información adicional separada por puntos y coma.

| **Campo**   | **Significado**                                              | **Descripción**                                                                                                                                                      | Programa |
| ----------- | ------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- |
| **DP**      | Depth                                                        | Profundidad de lectura aproximada de la variante                                                                                                                     |          |
| AS_SB_TABLE | Strand Bias Table                                            | Conteos (directa y reversa) para cada alelo (separados por una *pipe*)                                                                                               | Mutect2  |
| ECNT        | Expected alternate allele count                              | Calculation of the expected number of alternate alleles in a sample, based on the allele frequency in a population and the total number of reads at a given position |          |
| MBQ         | Mapping Quality of the Bases                                 | La calidad de las bases en una posición determinada del genoma.                                                                                                      |          |
| MMQ         | Median Mapping Quality (of the bases supporting each allele) | The median mapping quality of the reads that support each allele at a given position in the genome.                                                                  |          |
|   MPOS          |                Median Position                                           | Median position of the variant in the read                                                                                                                                                                      |          |

## FORMAT
