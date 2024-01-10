---
tags:
  - phd
  - write
---
Primero, voy a haceros un repaso global de los principales estudios y métodos que he encontrado, y después me gustaría que entre todos discutiéramos algunos puntos clave sobre la metodología que vamos a emplear para esta parte de epigenómica.

NOTA: Si algún estudio pone algo de correlación entre array y WGBS, ponerlo.
## Estudios detrás de la base

### Resistencia a fármacos

[[Methylation and drug resistance]]

### Enfermedad mínima residual

![[Metilación del ADN#^677262]]


## Marcadores encontrados en CRC

![[ctDNA methylation in CRC — biomarkers]]
## Principales metodologías y estudios

[[Revisión de estudios de metilación en ctDNA]]
### **WGBS**

El gold-standard en el análisis de metilación, porque permiten analizar citosinas únicas con gran precisión. Sin embargo, tienen como principales inconvenientes el daño al DNA y la baja cobertura.

En este punto, es importante destacar una importante diferencia de los estudios de WGBS con respecto al nuestro. ![[Research caveats#^162ba4]]
### cfMeDIP

### Arrays de metilación

## Aproximaciones bioinformáticas

Como ya os adelanté, el tratamiento bioinformático es muy similar al que se realiza para muestras de tejido. El pipeline lo estaría construyendo en base a [[Genomic benchmarking resources#SEQC2 Epigenomics (EpiQC)]].

Sin embargo, hay particularidades que se pueden utilizar para aumentar la señal del ctDNA frente al cfDNA. Os las comento un poco para que tengamos en mente las posibilidades de lo que podemos hacer.

### Detección de composición de tipos celulares

Existen diferentes programas para la deconvolución de datos epigenómicos en cuanto a la contribución de los tipos celulares que los conforman. 

De acuerdo con el principio de que los tumores provocan mayor liberación de cfDNA en el tejido subyacente, un estudio ha observado que se puede utilizar esta aproximación para identificar la presencia de cáncer en el ctDNA [@sunPlasmaDNATissue2015]. Estudio de 2015, baja cobertura (pero no la necesitan).

![[Pasted image 20231212103313.png]]

### Detección de firmas asociadas a cáncer

- CancerDetector
- CancerLocator

A parte, algunos otros artículos han seguido una estrategia por la que primero comparan datos epigenómicos de tejido (por ejemplo, disponibles en el TCGA) frente a sano para identificar las principales regiones de metilación diferencial específicas de CCR, y de ahí buscarlas en los datos de cfDNA.

## Puntos importantes

- Recapitular el tema de la cobertura y la calidad de las muestras. Me preocupa la cobertura que vayamos a obtener con WGBS. Tenemos que tener en cuenta que vamos a trabajar en contextos de enfermedad mínima residual, y que el tratamiento con bisulfito también daña las muestras.
- La RRBS tenía sentido porque nos permitía obtener más cobertura en zonas ricas en CpGs, que son lo que realmente nos importa. El genoma global tiene bajo contenido en CpGs, que se concentran en los promotores. WGBS tiene más sentido para estudiar niveles de metilación globales, más que estudiar citosinas individuales en este contexto.
- La aproximación que elijamos debe ser coste-efectiva. Si no, plantearnos una estrategia alternativa: MeDIP-seq, microarray.