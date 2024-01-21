---
tags:
  - phd
  - review
description: Esta revisión contiene información sobre todo el estado del arte del estudio de la metilación en el ctDNA.
---
The majority of ctDNA analysis focus on genomic alterations, which gives us only a partial picture of the tumor landscape [@bacaLiquidBiopsyEpigenomic2023]. Additional biomarkers are thus needed to improve ctDNA and liquid biopsy.

Aberrant methylation is implied in cancer pathogenesis. In fact, in most cancers, aberrant methylation is believed to be one of the early events in initiation. Epigenetic alterations are estable in ctDNA, and thus, can be a potential biomarker for liquid biopsy [@evilianidouDetectionRelevanceEpigenetic2021]. Similar to how mutations affect gene functions, methylation affects gene expression [@evilianidouDetectionRelevanceEpigenetic2021].

- **Early diagnosis.**
- **Detection of resistance mechanisms.** [[Methylation and drug resistance]].

Con respecto a las mutaciones, la metilación del ctDNA presenta características diferenciales clave:

- Ocurre a lo largo de regiones (islas CpG) vs. posiciones específicas
- Cambios universales, menor heterogeneidad intertumoral
- Mayor frecuencia general que mutaciones

Que han hecho que estos biomarcadores hayan suscitado un enorme interés en la detección temprana del cáncer y la detección de la EMR.


## Técnicas de detección de la metilación

### WGBS

Whole-Genome Bisulfite Sequencing (WGBS) is the most comprehensive and informative DNA methylation profiling technology. The **major advantage** is that it can detect the methylation state of every citosine, including low CpG density regions and non-CpG sites (CpA, CpT, and CpC) [@huangCellFreeDNAMethylation2019]. The **major disadvantage** is the low coverage and its **high cost**.

La incidencia de la metilación aberrante de ciertas islas CpG es relativamente elevada en las muestras tumorales, y por tanto puede ser más fácilmente identificada por técnicas de genoma completo como el WGBS, a diferencia de las mutaciones que no solo son raras, sino que se encuentran distribuidas en posiciones específicas.

Ventajas

- Tiene un gran potencial exploratorio, con una resolución de nucleótidos únicos.

Desventajas

- Produce degradación del ADN utilizado, lo que resulta en una menor cantidad de cfDNA utilizable para el análisis.

Studies: [[ctDNA methylation studies — methodologies and results#WGBS]].

### RRBS

### scWGBS y scRRBS

### Cell-free methylated DNA immunoprecipitation (cfMeDIPseq)

This approach allows for sensitive detection of cfDNA from small cfDNA quantities, and is cost-effective, more son than WGBS [@chenCellfreeDNAMethylome2022; @shenSensitiveTumourDetection2018].  

Studies: [[ctDNA methylation studies — methodologies and results#cfMeDIP]].

### Arrays de metilación

### Methylation-specific PCR

MSP is based on the use of two distinct methylation-specific primer sets for detecting the DNA of interest. The methylated primer will amplify bisulfite converted methylated DNA and untreated DNA, while the unmethylated primer is specific for bisulfite converted DNA in an unmethylated condition [@huangCellFreeDNAMethylation2019].

## Relevancia clínica en los procedimientos oncológicos
### Early diagnosis

Aberrant DNA methylation, in comparison with DNA mutations, has a series of specific advantages that make them more suitable for early cancer detection [@royDiagnosticPowerDNA2020]:

- They happen early in tumorigenesis, and can be tissue- and cancer-specific. 

- In regard to the tissue specificity, patterns of methylation in ctDNA can be used through deconvolution algorithms in order to determine the tissue of origin.

- DNA methylation patterns are widespread across the tumor tissue and across same tumor types. Somatic mutations, instead, are often limited to subpopulations or specific clones of tumor cells.

- DNA methylation is consistent across a large genomic region. This enables the use of multiple CpG sites for detection.

Por esta razón, la mayoría de los estudios que se han realizado han tenido enfoque en la detección temprana del cáncer.

A good review of current diagnostic tests can be found at: [@royDiagnosticPowerDNA2020]. Examples are:

- Epi proColon (see [[ctDNA methylation in CRC — biomarkers]]).
- Cologuard
- (...)
### Tumor-specific methylation detection

As I well know, the ctDNA concentration in cfDNA is generally low in cancer patients. The most promising approach to detect the specific ctDNA methylation patterns is to recover the tumor signal by using **deconvolution algorithms**. 

To identify ctDNA-specific methylation, one common strategy is to use the methylation profile of tumor-free peripheral blood mononuclear cells (PBMCs) as a negative control [@huangCellFreeDNAMethylation2019].  

Consult more studies in [@huangCellFreeDNAMethylation2019]. #check

### Detection of the tissue-of-origin

**The pattern of cfDNA methylation is consistent with their originated cells or tissues [@huangCellFreeDNAMethylation2019, @lianidouDetectionRelevanceEpigenetic2021].** Therefore, it can be used to infer the underlying cell type via **deconvolution methods**, and could be used to infer cancer subtypes (with implications in prognosis and treatment). #insight

### Detection of Minimal Residual Disease (MRD)

Hasta la fecha, la biopsia líquida centrada en la detección de [[Enfermedad mínima residual (EMR)]] se ha centrado en la detección de mutaciones somáticas. Sin embargo, el uso de la metilación aberrante del ctDNA como biomarcador para la detección de EMR es de enorme utilidad porque [@johnstonEpigeneticLiquidBiopsies2023a]:

- Los cambios epigenéticos son más frecuentes y universales que las mutaciones.
- El ctDNA mantiene muchos de estos cambios.

Epigenetic changes are more frequent and universal than genetic alterations in cancer, and ctDNA retains much of these changes, therefore making them suitable for MRD detection [@johnstonEpigeneticLiquidBiopsies2023a]. This key characteristic **holds great value for tumor-agnostic assays**; in fact, some of the most relevant commercial assays (LUNAR-1, etc.) incorporate methylation biomarkers.

A study in 51 metastatic cancer patients (INSPIRE trial, NCT02644369) showed that cancer specific methylation (CSM, using cf-MeDIPseq) predicted OS and PFS comparably to mutation concentration (MC) in ctDNA, both tumor-naïve approaches. Adding in the short fragment fraction (SFF) determination to the CSM, OS and PFS could be predicted better than the gold-standard MC. [Redirecting](https://doi.org/10.1016/j.annonc.2022.07.1744). ^31c690

Un estudio en [[cáncer de mama]] localizado identificó como la detección de cfDNA a través de los patrones de metilación del cfDNA con una sensibilidad del 80% y una especificidad del 97%.

### Respuesta a terapia y resistencia

## Clinical trials

Una buena revisión de los ensayos clínicos activos a abril de 2023 se encuentra en: [@johnstonEpigeneticLiquidBiopsies2023a].

- ECLIPSE (Guardant Health)
- MethylGene Tech
- GRAIL
- Clinical Genomics

## Estudios específicos

![[ctDNA methylation studies — methodologies and results]]