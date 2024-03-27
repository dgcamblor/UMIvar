---
aliases:
  - Bisulfite sequencing
---
La **secuenciación por [[Conversión por bisulfito]]** se considera el *gold standard* del estudio del metiloma.

## Whole-genome bisulfite sequencing (WGBS)

La secuenciación de genoma completo por conversión de bisulfito (WGBS) representa la mejor forma de obtener información del metiloma al completo a nivel de cada sitio CpG individual. Esta técnica consiste en la secuenciación de genoma completo por NGS de librerías que incorporan un paso de conversión por bisulfito.

En comparación con los [[Arrays de metilación]]:

- La ventaja es obvia: se obtiene información de muchos más sitios e incluso sitios no CpG (en el caso de que se necesiten estudiar). Se estima que puede llegar a cubrir aproximadamente el 95% de todas las citosinas en los genomas conocidos.
- La desventaja es también obvia: se requieren unos costes muy elevados para perfilar todo el genoma, especialmente a una cobertura que sea de utilidad (~30X). Se beneficia, por supuesto, del abaratamiento de los costes de secuenciación.

## Reduced-Representation Bisulfite Sequencing (RRBS)

La RRBS representa una variación con respecto al protocolo de WGBS que permite enriquecer las regiones analizadas en sitios más informativos para un análisis más coste-efectivo. En el protocolo más básico, durante la preparación de la librería se integra un paso de digestión con el enzima de restricción *Msp*I.

La enzima de restricción MspI corta el DNA en todos los sitios 5′-CCGG-3′ (independientemente del estado de metilación). Esta secuencia es muy frecuente en regiones del genoma ricas en CpG, como las islas CpG. El paso que permite el enriquecimiento es la selección de fragmentos de DNA en un rango de 40-220 pares de bases. Tras la digestión, los sitios ricos en CpG tienden a producir fragmentos de menor tamaño.

![[Pasted image 20240327195734.png]]