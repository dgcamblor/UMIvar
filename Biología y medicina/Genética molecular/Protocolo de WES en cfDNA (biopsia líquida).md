---
tags:
  - phd
---

## Extracción del DNA

Se toma la sangre del paciente. Dado que el cfDNA se encuentra en la fracción líquida de la sangre, se necesita realizar un protocolo de separación de la sangre:

- **Plasma:** fracción líquida de la sangre, con factores de coagulación (fibrinógeno), sin células (glóbulos blancos, plaquetas, glóbulos rojos). Se obtiene por centrifugación de un tubo de sangre en presencia de anticoagulante.
- **Suero:** fracción líquida de la sangre, sin factores de coagulación y sin células. Se obtiene por centrifugación de un tubo de sangre que se ha dejado coagular de forma natural.

El plasma es preferido para la extracción de cfDNA [@pittella-silvaPlasmaSerumWhich2020].

El ADN debe ser extraído del plasma utilizando un kit específico de extracción de cfDNA (también para otro DNA/RNA en sangre):

- [QIAamp Circulating Nucleic Acid kit](https://www.qiagen.com/es-us/products/discovery-and-translational-research/dna-rna-purification/dna-purification/cell-free-dna/qiaamp-circulating-nucleic-acid-kit)
- Etc.

## Cuantificación del DNA

Es necesario comprobar que la cantidad y el estado del DNA extraído sea la óptima. Para ello, se necesita realiza un ensayo basado en electroforesis, como el [Cell-free DNA ScreenTape Analysis](https://www.agilent.com/en/product/automated-electrophoresis/tapestation-systems/tapestation-dna-screentape-reagents/cell-free-dna-screentape-analysis-294915) para sistemas Agilent TapeStation.

De esta forma, se obtiene un perfil como el siguiente:

![[Profile_figure 2_2019-07-24.3_B1.webp]]

Donde las diferentes bandas se corresponden con:

- Primera banda (o la que corresponda según su tamaño). El ladder de DNA que se utiliza como control positivo, sirviendo como referencia para determinar el tamaño (en bp) de las bandas subsiguientes.
- Banda ~170 bp. Se corresponde con los fragmentos de [[cfDNA y ctDNA]] mononucleosomales (esto es, asociados a una única histona).
- Bandas subsiguientes. Se corresponden con fragmentos de DNA, menos abundantes, asociados a multímeros nucleosomales. 

Esto permite obtener un porcentaje de contenido de cfDNA de la muestra extraída. Se fija un porcentaje umbral, por debajo del cuál una muestra no es aceptable (por ejemplo, 70%).

Todo el cfDNA extraído se almacena en refrigeración.

## Preparación de las librerías para NGS-WES

Se extrae una alícuota del cfDNA extraído para la preparación de las librerías. Suelen utilizarse 10-40 ng de cfDNA.

Los pasos de la preparación de librerías son los habituales.

![[Preparación de librerías para NGS]]

Posteriormente, la librería se enriquece en exoma utilizando un kit de captura, como [KAPA HyperExome kit](https://sequencing.roche.com/global/en/products/group/kapa-hyperexome.html).

La preparación de librerías para NGS-WES requiere de un control de calidad tanto en la librería **pre-captura** como en **post-captura**, lo que puede realizarse mediante ensayo basado en electroforesis (como [Cell-free DNA ScreenTape Analysis](https://www.agilent.com/en/product/automated-electrophoresis/tapestation-systems/tapestation-dna-screentape-reagents/cell-free-dna-screentape-analysis-294915)).
## Secuenciación

Los principales equipos de secuenciación disponibles en el mercados son los de Illumina. Una comparación entre equipos se encuentra disponible [aquí](https://emea.illumina.com/systems/sequencing-platforms/comparison-tool_msm_moved_msm_moved.html).

- HiSeq 3000
- NovaSeq6000

Al igual que en otro tipo de muestras, la [[Cobertura]] objetivo debe ir acorde con la frecuencia a la que se encuentran las alteraciones moleculares que se deben encontrar.