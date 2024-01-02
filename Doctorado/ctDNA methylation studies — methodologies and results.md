---
tags:
  - phd
  - review
description: Specific studies that have used WGBS. A comment on their approach. Writing down the methodologies they use.
---

The review [@johnstonEpigeneticLiquidBiopsies2023] contains several clinical trials in development.

## WGBS

- Whole-genome bisulfite sequencing of cell-free DNA identifies signature associated with metastatic breast cancer. [@legendreWholegenomeBisulfiteSequencing2015].
	- Differential methylation analysis between healthy individuals (H), disease-free survivors (DFS) and metastatic breast cancer (MBC) patients.
	- cfDNA methylomes of MBC patients are different from disease-free survival and healthy controls. They define 21 DNA hypermethylation hotspots.

- Whole-genome circulating tumor DNA methylation landscape reveals sensitive biomarkers of breast cancer. [@haiWholeGenomeCirculating2022]. 
	- They propose an improved way to extract ctDNA (in order to reduce background contamination). Patients of breast cancer (n = 74) vs. age-matched healthy controls (n = 7). 26 CpG sites that could serve as potential biomarkers of diagnosis ^af677d
	- **Pipeline:** [GitHub - zhq921/cfWGBS-bioinfo-pip: A comprehensive and rigid computational framework to identify reliable cfDNA methylation biomarkers from cfWGBS data](https://github.com/zhq921/cfWGBS-bioinfo-pip)
	- **Data:** [Browse - GSA - CNCB-NGDC](https://ngdc.cncb.ac.cn/gsa/browse/CRA001142)

- Noninvasive detection of cancer-associated genomewide hypomethylation and copy number aberrations by plasma DNA bisulfite sequencing. [@chanNoninvasiveDetectionCancerassociated2013]. 
	- Referenced by [@shenSensitiveTumourDetection2018]. 
	- Analyzes genome wide (hypo)methylation as a biomarker for hepatocelular carcinoma. Not detecting specifical biomarkers. Bins of 1 Mb were considered.

- [@sunPlasmaDNATissue2015]. 
	- Referenced by [@shenSensitiveTumourDetection2018]. Deconvolution of genome-wide bisulfite sequencing of plasma DNA. They study the contribution of each tissue to cfDNA in several conditions, including HCC patients. The contribution of liver tissue to cfDNA is higher in HCC patients than in controls. 

- **Hypomethylation in HBV integration regions aids non-invasive surveillance to hepatocellular carcinoma by low-pass genome-wide bisulfite sequencing.** [@zhangHypomethylationHBVIntegration2020]. 
	- They propose the strategy of **low-pass WGBS** to detect methylation changes in circulating cell-free DNA (cfDNA) from patients with liver diseases and hepatocellular carcinoma (HCC). The challenge of cfDNA is to specifically detect the ctDNA patterns, as it is a minor fraction of cfDNA, specially in early-stage cancer and MRD settings, something which benefits from deep sequencing. However, this study proposes a reduced depth sequencing of high sample sizes for cohort studies.
	- Demuestran que la estrategia de low-pass WGBS es eficaz 

- [@luCellfreeDNACfDNA2019]. NOTE: Meeting abstract. No paper found.
	- Full text was not available.
	- Mediante WGBS, identifican 729 DMRs entre pacientes que recurren y que no recurren. Hacen clustering no supervisado usando dichas DMRs para distinguir entre los grupos de recurencia.

- Whole-genome bisulfite sequencing analysis of circulating tumour DNA for the detection and molecular classification of cancer [@gaoWholeGenomeBisulfite2022].
	- Describen un método de análisis de ctDNA a partir de pequeñas muestras de plasma (200 uL). 

- Methylation of NBPF1 as a novel marker for the detection of plasma cell-free DNA of breast cancer patients. [@liMethylationNBPF1Novel2018].
	- Comparan 25 pacientes de cáncer de mama, 25 con enfermedad de mama benigna y 25 controles sanas.
	- Encuentran metilación diferencial en NBPF1, y confirman su potencial diagnóstico de este gen para el cáncer de mama.

- Targeted methylation sequencing of plasma cell-free DNA for cancer detection and classification. [@liuTargetedMethylationSequencing2018].
	- Identificación de 9223 sitios CpG hipermetilados en pan-cáncer según datos de TCGA. Los usaron para la correlación con diferentes tipos de cáncer. El score de metilación correlaciona con el resultado del tratamiento.
	- WGBS seguida de enriquecimiento en dianas.

> [!NOTE] Summary
> WGBS has been used specifically for cancer detection.

## cfMeDIP

- [@chenCellfreeDNAMethylome2022]. 
	- Analysis of the whole cfDNA methylome in 60 plasma samples from localized prostate cancer vs. 175 plasma samples from metastatic prostate cancer. Uses cfMeDIP, based on the advances made by Shen et al [@shenSensitiveTumourDetection2018].
- [@shenSensitiveTumourDetection2018].

## Methylation arrays

- [Identification of DNA methylation biomarkers with potential to predict response to neoadjuvant chemotherapy in triple-negative breast cancer](https://doi.org/10.1186/s13148-021-01210-6). 
	- Prediction of response to neoadjuvant therapy in triple-negative breast cancer patients, prior to resection. Identification of differentially methylated probes (DMPs) between complete, partial and non-responders.

- [@luoCirculatingTumorDNA2020b]. 
	- NOTE: Editorial Concern. Prospective study of screening for high risk CRC patients. They elaborate a prognostic prediction machine learning model to predict risk for CRC. They detected [[Detection of CRC-specific methylation profiles | CRC-specific methylation markers]].

- Liquid biopsy epigenomic profiling for cancer subtyping. [@bacaLiquidBiopsyEpigenomic2023]

- A new approach to epigenome-wide discovery of non-invasive methylation biomarkers for colorectal cancer screening in circulating cell-free DNA using pooled samples. [A new approach to epigenome-wide discovery of non-invasive methylation biomarkers for colorectal cancer screening in circulating cell-free DNA using pooled samples | Clinical Epigenetics | Full Text](https://doi.org/10.1186/s13148-018-0487-y).

## ddPCR

- [@boeckxMutationMethylationAnalysis2018]. They follow 24 mCRC patients using ddPCR for selected mutations and methylation of the *NPY* gene. This gene has been found to be hypermethylated in CRC tissue.