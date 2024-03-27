---
aliases:
  - Methylated DNA immunoprecipitation
  - MeDIP
---
La inmunoprecipitación del ADN metilado (MeDIP) consiste en la fragmentación del ADN (el tamaño de los fragmentos marca la resolución) y el posterior aislamiento de los fragmentos que contienen ADN metilado via anticuerpos específicos anti 5-metilcitosina o 5-hidroximetilcitosina.

![[Pasted image 20240327201733.png]]

Esta técnica cubre CpGs (y no CpGs) distribuidos a lo largo del genoma, en todo tipo de regiones (densas, no densas, repetitivas). En comparación con las técnicas basadas en [[Conversión por bisulfito]], el uso de anticuerpos específicos permite seleccionar la modificación concreta que se quiere estudiar: 5-metilcitosina o 5-hidroximetilcitosina. A cambio, la resolución es mucho menor (~150 bp, la de los fragmentos inmunoprecipitados.)

![[Pasted image 20240327201038.png]]

## MeDIP-seq

Cuando la técnica MeDIP se complementa con NGS, se conoce como MeDIP-seq. Mediante la secuenciación de los fragmentos de ADN inmunoprecipitados y s u mapeo al genoma, la cobertura de las regiones se utiliza para estimar el nivel de metilación de cada una de ellas.

- Cuando una región tiene mayor cobertura, esto es que se han inmunoprecipitado más fragmentos, lo que a su vez depende de que tenga elevados niveles de metilación.
- Y viceversa.

![[Pasted image 20240327202318.png]]