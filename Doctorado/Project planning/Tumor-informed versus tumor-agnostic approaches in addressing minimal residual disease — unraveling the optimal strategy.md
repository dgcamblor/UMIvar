---
tags:
  - paper
  - phd
---
TP -> ctDNA+ and recurring
FP -> ctDNA+ and not recurring
TN -> ctDNA- and not recurring
FN -> ctDNA- and recurring

---

$$
Sensibilidad = TP/(TP+FN)
$$
$$
Especificidad = TN/(TN+FP)
$$
$$
VPP = TP/(TP+FP)
$$
$$
VPN = TN/(TN+FN)
$$

---

$$
Assay sensitivity = Recurring\ ctDNA+/\ Total\ recurring
$$
$$
Assay specificity = Non-recurring\ ctDNA-/\ Total\ non-recurring
$$
$$
VPP = Recurring\ ctDNA+/\ Total\ ctDNA+
$$
$$
VPN = Non-recurring\ ctDNA-/\ Total\ ctDNA-
$$

---

Quizás deberíamos de desarrollar este párrafo un poco más porque veo que el 91% de sensibilidad ese es un poco muy polémico. En el estudio de Parikh, reportan:  
  
- Landmark: 55.6% sensibilidad, 100% especificidad.  
- Longitudinal: 69% sensibilidad, 100% especificidad.  
- Surveillance (metodología de Reinert et al., considerando muestras sacadas en un rango de 4 meses de la recurrencia): 91% sin reporte de especificidad.

Natera y Guardant Health entraron en disputas legales sobre este punto, de si es misleading etc., dejo aquí una contestación de Natera: [https://www.natera.com/company/news/natera-responds-to-meritless-lawsuit-and-files-false-advertising-claim-against-guardant-for-misleading-oncologists-2/](https://www.google.com/url?q=https://www.natera.com/company/news/natera-responds-to-meritless-lawsuit-and-files-false-advertising-claim-against-guardant-for-misleading-oncologists-2/&sa=D&source=docs&ust=1711445708106011&usg=AOvVaw2fxanPE3uYK-80L6UcUb_s)
Podemos poner ambos valores de sensibilidad en la tabla (69% y 91%) y poner este último con asterisco por las particularidades con las que se ha calculado.

Regardless of the technique, ctDNA assays for the detection of MRD will require high specificity to escalate adjuvant therapy in patients who would be categorized as low risk according to clinicopathologic factors and would ultimately recur without additional therapy, and high sensitivity to allow safe de-escalation of therapy in those patients who are traditionally treated with intensive chemotherapy and avoidable adverse events. [@jacomeMinimalResidualDisease2023]. 
> Un poco a colación de lo que dicen Reinert 2024

>In a prospective study with 112 mCRC patients who had undergone liver resection with curative intent, ctDNA positivity by a tumor-informed assay was the most significant prognostic factor for disease-free survival (HR: 5.7, 95% CI 3.3–10.0) [52]. MRD was detected in 61 (54%) patients, of which 59 (97%) presented recurrence at the time of data cutoff.

## II. Tumor-agnostic strategies

### Methylation

[Efficient detection and post-surgical monitoring of colon cancer with a multi-marker DNA methylation liquid biopsy | PNAS](https://www.pnas.org/doi/full/10.1073/pnas.2017421118)

The use of techniques that detect epigenomic signatures assumes that hypermethylation of promoter regions of some tumor suppressor genes is a common finding in CRC. Hypermethylation of WNT-inhibitor-factor-1 (WIF1) and neuropeptide Y (NPY) has been validated as a marker of both advanced and localized colon cancer [35,36].

[Response prediction by mutation- or methylation-specific detection of ctDNA dynamics in pretreated metastatic colorectal cancer - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10541738/)

## III. Tumor-informed strategies

Importance of ctDNA surveillance, % of patients that still relapse after postoperative (1 month) ctDNA-negative status: 9.1% [@chenPostoperativeCirculatingTumor2021], 11.9% [@reinertAnalysisPlasmaCellFree2019].

[Comparing single‐target and multitarget approaches for postoperative circulating tumour DNA detection in stage II–III colorectal cancer patients - Henriksen - 2022 - Molecular Oncology - Wiley Online Library](https://febs.onlinelibrary.wiley.com/doi/10.1002/1878-0261.13294) -> Comparison between ST and MT strategies.
[Study Details | | ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT03748680)

>Chan et al. “After the exclusion of CH mutations, 37% (14/38) of the evaluated patients harbor at least one tumor-derived mutation from plasma cfDNA for longitudinal monitoring.” “In contrast to the tumor-agnostic approach, up to 84% (32/38) of the patients harbor at least one mutation from the tumor tissue for subsequent plasma cfDNA monitoring (Figure 3A).”  Interestingly, ctDNA was detected in all three metastatic patients from our study cohort using both the tumor-informed and tumor-agnostic approach, highlighting that the impact of assay sensitivity is more prominent in patients with localized CRC.

>Another direct comparison although at a smaller scale came from a prospective study by Chan et al....

### ddPCR 

### SafeSeqS approach: the DYNAMIC studies

> [Program Guide – ASCO Meeting Program Guide](https://meetings.asco.org/abstracts-presentations/228859) 
> [Studies Confirm the Utility of ctDNA in Guiding Adjuvant Chemotherapy in CRC | ASCO Daily News](https://dailynews.ascopubs.org/do/studies-confirm-utility-ctdna-guiding-adjuvant-chemotherapy-crc)

>

### WES in tissue characterization: Signatera

- MEDOCC-CrEATE trial
	- Protocol: [Circulating tumor DNA guided adjuvant chemotherapy in stage II colon cancer (MEDOCC-CrEATE): study protocol for a trial within a cohort study | BMC Cancer | Full Text](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07252-y)
	- Currently enrolling: search for more

> Each target can achieve an average depth of 100,000x with ultra‑deep sequencing, and this approach has exhibited a sensitivity >98% at 0.01‑0.02% ctDNA concentrations (60). [@chenDetectingLiquidRemnants2023] The predefined threshold is based on Natera’s proprietary variant calling method wherein detecting at least 2 out of 16 variants ensures the optimal analytical performance of the assay with > 95% sensitivity at 0.01% mean variant allele frequency and with 99.7% specificity [@kotaniMolecularResidualDisease2023].


## IV. Comparison

- Tumor-informed assays demand more time than tumor-naïve approaches for personalization.
- As the tumor-informed approach employs a tissue sample for variant analysis, this enables the filtering of CHIP mutations found in cfDNA (decreasing false positive rates).
- Tumor-informed assays are reported to have greater sensitivity than tumor-naive panels, even when they use large gene panels [@kasiImpactCirculatingTumor2022]. A tumor-naïve panel may not cover variants found in some patients.
- Currently, the development of tumorinformed assays is challenging in the clinical setting given the high cost, the mandatory requirement for tumor and germline samples, and the time required for patient-specific assay development [@santonjaComparisonTumorinformedTumornaive2023].

The body of research in CRC:

- Circulating Tumor DNA As A MRD Assessment And Recurrence Risk In Patients Undergoing Curative Intent Resection With Or Without Adjuvant Chemotherapy In Colorectal Cancer: A Meta-analysis [@chidharlaCirculatingTumorDNA2022].
	- Meta-analysis of studies detecting ctDNA as a proxy for MRD in Stage I-IV CRC patients after curative-intent surgery.
	- They perform a subgroup analysis for tumor-informed vs. tumor-agnostic studies, favoring the former over the later.

- Tumor-informed or tumor-agnostic circulating tumor DNA as a biomarker for risk of recurrence in resected colorectal cancer patients [@chanTumorinformedTumoragnosticCirculating2023].
	- They use both a tumor-agnostic and a tumor-informed approach to determine their sensitivity.
	- Utilizing the tumor-agnostic approach reduced the recurrence detection sensitivity to 67%.

Studies in other types of cancer:

- Comparison of tumor-informed and tumor-naïve sequencing assays for ctDNA detection in breast cancer [@santonjaComparisonTumorinformedTumornaive2023]. 
	- Uses whole-genome rather than whole-exome for somatic mutation identification in the tumor-informed approach.

## V. Emerging Technologies and Future Directions

- [x] Añadir metodología de MAESTRO ✅ 2024-03-27

### 1. Breakthroughs in ctDNA genomics

> [FoundationOne Tracker | Foundation Medicine](https://www.foundationmedicine.com/test/foundationone-tracker) [Study Details | | ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT04259944)



[CODEC]

[MAESTRO]

 lowering of sequencing costs, the emergence of new technologies for low-frequency variant enrichment (e.g., MAESTRO [36]) and background noise error correction and the developments in the field of bioinformatics.

### 2. Epigenomics and beyond

> Incorporation of other molecular markers cfRNA, miRNA...

> ¿Añadir falsos positivos?

> WGBS/cfMEDIP. In the whole-genome setting, the gold standard for methylation assessment, Whole Genome Bisulfite Sequencing

- [x] Discutir el tema de los falsos positivos en Guardant Reveal. ✅ 2024-03-27

### 3. Potential new biomarkers in liquid biopsy

### 4. Advanced Imaging Techniques

> The application of AI in radiomics.
> However, ctDNA still takes 


### 4. (Bioinformatic methods and) AI

> Beyond radiomics, AI also has interesting applications...
> (CASTLE,19 Poisson,19,20 ALPACA,19,21 and Dynamic LOB19) [@henriksenUnravelingPotentialClinical2024]

## Challenges and future directions

> False negatives and positives, landmark vs sequential...

> Importance of the algorithm of ctDNA calling. Which is the cutoff for saying a patient is ctDNA positive?
