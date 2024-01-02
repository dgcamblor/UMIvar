---
tags:
  - phd
---


El CCRm constituye un importante problema clínico dentro del cáncer colorrectal, donde la mediana (Me) de supervivencia de los pacientes es de 36 meses.

- 15-30% pacientes son diagnosticados en etapa avanzada.
- 20-50% desarrolla metástasis tras enfermedad localizada.

A pesar de ello, se ha habido pocos avances en tratamiento que tengan un impacto en el pronóstico.

La biopsia líquida, que como sabéis ha demostrado enormes ventajas en la caracterización tumoral con objetivos terapéuticos en CCR localizado, sería una importante herramienta que nos podría ayudar a identificar oportunidades de tratamiento en CCRm.

La biopsia líquida permite, a nivel genómico, obtener un panorama más general de las mutaciones que tiene un tumor, atendiendo a las diferentes subpoblaciones de células tumorales, frente a la biopsia del tejido tumoral, en la que solo se caracteriza la sección biopsiada. Sin embargo, estudiar el panorama mutacional es quedarnos con una visión parcial de la dinámica tumoral en la metástasis, puesto que en ella participan muchos otros procesos.

Es en este punto donde entramos a considerar la epigenética, que comprende una serie de mecanismos dinámicos que regulan la expresión de los genes. En concreto, la metilación de citocinas en su carbono es una marca que ha sido altamente relacionada con todos los mecanismos tumorales desde la tumorigénesis hasta la progresión a enfermedad avanzada. Esta marca se encuentra de forma estable en el ADN tumoral circulante y, por lo tanto, puede ser utilizada como parte del análisis de la biopsia líquida. Sin embargo, aunque se han hecho avances en su uso, su uso en la práctica clínica requiere de más estudios que demuestren su beneficio.

Con todo esto, la hipótesis de nuestro proyecto es que la integración de los datos genómicos mutaciones con los datos epigenómicos permitirá obtener una perspectiva sobre la dinámica tumoral en la metástasis y revelar potenciales dianas accionables.

Objetivos.

Para nuestro estudio, trabajaremos con una cohorte de 36 pacientes de CCRm, para los que se hará un análisis de secuenciación de exoma completo (WES) y secuenciación de bisulfito de genoma completo (WGBseq) en diagnóstico y progresión tumoral. De esta forma, se localizarán las principales mutaciones y sitios CpG implicados. Estas alteraciones moleculares se seguirán en la respuesta al tratamiento mediante secuenciación de amplicones y PCR específica de metilación.

Por otra parte, a partir de estos resultados se hará el ensayo de nuevos fármacos en organoides.

Respecto a la parte bioinformática, para la parte genómica emplearé el *pipeline* de rutina implementando el cambio que hemos comentado (fgbio) a partir de los resultados que he obtenido en mi TFM, ya que he visto que aumenta tanto la sensibilidad como la precisión de la llamada de variantes.

Para el análisis de metilación del ctDNA se hará un *pipeline* desde 0, incorporando herramientas clásicas para el análisis de este tipo de datos como alineadores de secuencias convertidas por bisulfito, *software* de *methylation calling*, etc. 

## Insights

Los estudios del ctDNA se han centrado en la búsqueda de alteraciones genómicas, pero esto ha limitado la habilidad para detectar variables clínicas de enorme importancia como lo son la expresión de los genes (Baca et al., 2023).

Es por eso que es tan importante incorporar la dimensión epigenética en este proyecto.

(Una mutación en un gen no tiene implicaciones si este ha sido silenciado epigenéticamente.)