---
tags:
  - phd
---
Tumor-informed and tumor-agnostic approaches are two different methods for detecting minimal residual disease (MRD) using circulating tumor DNA (ctDNA) as a biomarker. Tumor-informed assays use information from a patient's tumor sample to design a personalized test, while tumor-agnostic assays do not require prior tumor genomic knowledge and are more generic. Both approaches have their own advantages and limitations.

>A hot topic at the [American Society of Clinical Oncology (ASCO) 2023](https://conferences.asco.org/) in Chicago was approaches to liquid biopsy and how to best utilize them in clinical practice. One of the main debates is a liquid only (tumor naïve) or tissue plus liquid (tumor-informed) approach to testing.
>
>The **tumor-informed** approach involves using prior knowledge of the patient’s tumor tissue molecular characteristics to guide the analysis of liquid biopsy samples. Ideally, this tumor informed approach would be based off a baseline knowledge from a next-gen sequencing assay consisting of a targeted larger gene panel (>50 genes) or **whole genome sequencing** (WGS) that outline most common mutations in the sample. By incorporating this prior knowledge, the liquid biopsy assay can target specific genetic alterations associated with the tumor and provide more sensitive and specific results. However, due to the heterogeneity of tumor tissue and the clonal variation of tumor tissue during the treatment process, the initial specific tumor tissue genome cannot represent the genomic variation of the whole tumor tissue and the new tumor clonal variation.
>
>On the other hand, **tumor-agnostic assays**, which are independent of prior tumor genomic knowledge, may offer comparable sensitivity to tumor-informed assays and have the advantage of a more rapid turnaround time and reduced cost.
>
>Clinical studies have demonstrated increased sensitivity of tumor-informed MRD testing, especially for heterogeneous tumors.

## Outline final

I. Introduction
	1. Significance 

## Outline

**I. Introduction** 

1. Background 
2. Significance of MRD in Oncology 
3. Evolving Landscape of Tumor-Informed and Tumor-Agnostic Approaches

>Introduction to MRD, ctDNA and its potential usage for MRD detection, overview of ctDNA detection strategies (?).
>Detection of ctDNA can be particularly challenging in settings of low tumoral burden, such as early disease and minimal residual disease monitoring

**II. Tumor-Informed Strategies** 

**III. Tumor-Agnostic Approaches --> EC en marcha a nivel mundial empleando ambos enfoques**

**IV. Comparative Analysis ---> EC publicados (PEGASUS / LUNAR1, GALAXY / SIGNATERA, DYNAMIC / SAFESEQ)**

1. Strengths and Weaknesses of Tumor-Informed Approaches 
2. Strengths and Weaknesses of Tumor-Agnostic Approaches 
3. Integration of Complementary Strategies

>Incluir estudios específicos que resalten las diferencias entre una aproximación y otra.

**V. Emerging Technologies and Innovations** 

1. Next-Generation Sequencing. 
2. Liquid Biopsy Technologies 
3. Advanced Imaging Techniques

> Contar Next-Generation Sequencing en el primer apartado.
> Liquid Biopsy Technologies -> CTCs for MRD assessment?

**VI. Future Directions --> WGS, Epigenomics, fragmentomics**

1. Advancements in Precision Medicine 
2. Potential Therapeutic Implications 
3. Role in Personalized Treatment Paradigms

Tumor-naïve assays can't compete in terms of detection sensitivity with tumor-informed assays. However, the use of additional biomarkers could prove useful to enhance their capabilities. The use of ctDNA methylation could prove useful.

**VII. Conclusion** 

1. Recapitulation of Key Findings 
2. Recommendations for Future Research

## Cuestiones

- Parece que está claro que los assays tumor-informed presentan resultados más sensibles (porque están centrados en las mutaciones específicas del tumor concreto y no en regiones que pueden o no ser relevantes) y específicos (no se identifican falsos positivos en regiones no relevantes para el tumor). 
	- Esta observación está fundamentalmente basada en estudios individuales, pero hay pocos estudios que hayan hecho una comparación directa entre ellos [@santonjaComparisonTumorinformedTumornaive2023].
	- Sería bueno recopilar todos los estudios posibles (IV. Comparative Analysis) que hayan comparado objetivamente estas dos estrategias (ver [[Tumor-informed versus tumor-agnostic approaches in addressing minimal residual disease — unraveling the optimal strategy#Studies comparing both approaches]]). P. ej.: el metaanálisis que encontré.

- ¿Qué debería ir en el apartado de Next-Generation Sequencing, **V. Emerging Technologies and Innovations**? Muchos assays como Signatera ya están incorporando la NGS...
	- Quizás incorporar un apartado previo en **I. Introduction** sobre la descripción de las técnicas de detección de ctDNA actuales: similar a lo que se encuentra en [@zhuPlasmaAssayCellFree2023].
		- Métodos basados en PCR
		- Métodos basados en NGS
	- ¿Hace referencia a que NGS sea emergente? ¿O contar algo dentro de la NGS que sea emergente? Los UMIs ya llevan unos años... El propio DYNAMIC usa SafeSEQ

- En epigenética se puede hablar un poco de las particularidades y hasta ensayos clínicos que siguen una u otra aproximación.
	- Tengo una referencia en cáncer de pulmón de un estudio que comparó directamente las dos estrategias usando puramente marcadores epigenéticos
	
## Main advantages and disadvantages

- Tumor-informed assays demand more time than tumor-naïve approaches for personalization.
- As the tumor-informed approach employs a tissue sample for variant analysis, this enables the filtering of CHIP mutations found in cfDNA (decreasing false positive rates).
- Tumor-informed assays are reported to have greater sensitivity than tumor-naive panels, even when they use large gene panels [@kasiImpactCirculatingTumor2022]. A tumor-naïve panel may not cover variants found in some patients.
- Currently, the development of tumorinformed assays is challenging in the clinical setting given the high cost, the mandatory requirement for tumor and germline samples, and the time required for patient-specific assay development [@santonjaComparisonTumorinformedTumornaive2023].

## Studies comparing both approaches

%%TODO
- [ ] Study the nature of these studies. Use several headings as needed
- [ ] Add studies in epigenetic data
%%

### Genomics

- Circulating Tumor DNA As A MRD Assessment And Recurrence Risk In Patients Undergoing Curative Intent Resection With Or Without Adjuvant Chemotherapy In Colorectal Cancer: A Meta-analysis [@chidharlaCirculatingTumorDNA2022].
	- Meta-analysis of studies detecting ctDNA as a proxy for MRD in Stage I-IV CRC patients after curative-intent surgery.
	- They perform a subgroup analysis for tumor-informed vs. tumor-agnostic studies, favoring the former over the later.

- Comparison of tumor-informed and tumor-naïve sequencing assays for ctDNA detection in breast cancer [@santonjaComparisonTumorinformedTumornaive2023]. 
	- Uses whole-genome rather than whole-exome for somatic mutation identification in the tumor-informed approach.

- 

### Epigenomics

## Clinical trials

![[MRD ctDNA clinical trials]]

## Assays

![[MRD ctDNA assays]]
## References

- [Τumor-informed Vs Tumor-naïve | Genomic Testing Cooperative](https://genomictestingcooperative.com/tumor-informed-vs-tumor-naive/)