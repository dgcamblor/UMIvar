---
tags:
  - paper
  - phd
---
> In recently presented results, however, the [[COBRA trial]] was not able to meet its primary endpoint, and was promptly stopped due to futility in achieving statistical significance [@morrisPhaseIIResults2024].
## II. Tumor-agnostic strategies

> —also termed tumor-uninformed or tumor-naïve—

### Methylation

> Initial studies have explored the potential of single biomarker genes for methylation. (Then, GUARDANT)

[Response prediction by mutation- or methylation-specific detection of ctDNA dynamics in pretreated metastatic colorectal cancer - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10541738/)

## III. Tumor-informed strategies

- MEDOCC-CrEATE trial
	- Protocol: [Circulating tumor DNA guided adjuvant chemotherapy in stage II colon cancer (MEDOCC-CrEATE): study protocol for a trial within a cohort study | BMC Cancer | Full Text](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07252-y)
	- Currently enrolling: search for more

### DYNAMIC trial (SafeSeqS approach)

> Añadir sobre la capacidad de detección.
> El primero de TIE van solo a por una y el segundo por al menos una??? revisar

> Tumor-informed approaches have proven to be the most adopted 


A growing body of research and clinical trial results has revealed the usefulness of tumor-informed approaches to MRD detection. 

One of the earliest key studies that supported this approach was led by Tie et al., and was performed in a prospective cohort of 230 CRC stage II patients (trial ACTRN12612000326897), where ctDNA detection after resection showed higher risk of recurrence as determined by a CT scan. Initially, the variant with the highest allele frequency in FFPE tumor tissue was determined for each patient using targeted multiplex PCR aiming for 15 genes relevant to CRC cancer. Afterwards, those variants were tracked using a personalized Safe-SeqS assay, an UMI-based approach that enables deep sequencing with enhanced error correction [@kindeDetectionQuantificationRare2011]. Patients were then classified as ctDNA-positive or -negative, depending on the finding of said mutations in plasma. The results showed that postoperative ctDNA-positive status was related to higher risk of recurrence among patients not treated with adjuvant chemotherapy (HR = 18; 95%CI = 7.9-40) [@tieCirculatingTumorDNA2016].

> The growing evidence for the clinical utility of ctDNA in the postoperative setting encouraged the same group to conduct a series of related clinical trials aimed at evidencing the benefit of ctDNA determination to guide the adjuvant therapy. These prospective studies —DYNAMIC, DYNAMIC-III and DYNAMIC-Rectal— randomized patients with CRC stage II, stage III and locally advanced rectal cancer after receiving surgery into two groups: patients receiving treatment according to their ctDNA status and patients receiving conventional treatment. To determine the ctDNA presence, the studies adopted the aforementioned Safe-SeqS tumor-informed approach.


The DYNAMIC trial was the first to reveal insightful results 


The same tumor-informed approach was adopted in the DYNAMIC trial (ACTRN12615000381583), which further studied the possibility of using the postoperative ctDNA status as a basis for making treatment decisions for patients with CRC stage II. By random assignment, patients were assigned to an standard management group (n = 153) or a ctDNA-guided management group (n = 302). In this approximation, patients that were ctDNA-positive at 4 or 7 weeks after surgery were recommended to receive adjuvant oxaliplatin- or fluoropyrimidine-based chemotherapy. Conversely, patients testing negative for ctDNA did not receive adjuvant chemotherapy. The RFS for the ctDNA-guided group was comparable to that of the standard management group (93.5% vs. 92.4%), indicating the suitability of chemotherapy reduction in ctDNA negative patients [@tieCirculatingTumorDNA2022].

> DYNAMIC-II es el que se ha comentado anteriormente. Está también el DYNAMIC-III (ACTRN12617001566325) en pacientes de CCR III, y el DYNAMIC-RECTAL (LARC?).

### Other studies

Several other studies have explored the potential of assessing tumor mutations to inform the tracking of MRD via ctDNA.

![[Tumor-informed versus tumor-agnostic approaches in addressing minimal residual disease — unraveling the optimal strategy (planning)#Tumor-informed]]

A growing body of research has 

Studies using non-commercial tumor-informed assays have also rendered insightful results.



Other studies that have since relied on custom



[Studies Confirm the Utility of ctDNA in Guiding Adjuvant Chemotherapy in CRC | ASCO Daily News](https://dailynews.ascopubs.org/do/studies-confirm-utility-ctdna-guiding-adjuvant-chemotherapy-crc)
> Another similar technology is AlphaLiquid Detect [Personalised circulating tumour DNA assay with large-scale mutation coverage for sensitive minimal residual disease detection in colorectal cancer | British Journal of Cancer](https://www.nature.com/articles/s41416-023-02300-3).

[FoundationOne Tracker | Foundation Medicine](https://www.foundationmedicine.com/test/foundationone-tracker) [Study Details | | ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT04259944)

## IV. Comparison

### Strengths and weaknesses

- [ ] Añadir citas

As shown by the literature, tumor-informed approaches can detect MRD with higher sensitivity and specificity by looking for specific known mutations rather than relying on a panel of genes that may or may not be of enough relevance for the patient's tumor. This method also ensures that germline variants and variants related to clonal hematopoiesis of indeterminate potential are directly filtered out of the analysis [@steensmaClonalHematopoiesisIndeterminate2015]. Tumor-informed assays are, however, inherently limited by their reliance on tumor tissue. The most obvious examples in which this approximation is unfeasible are the cases where the tumor tissue is not accessible or the available sample is inadequate for genomic characterization, such as after extensive neoadjuvant treatment. Whenever the tumor tissue is available, the tissue biopsy and subsequent sequencing for the creation of a custom assay are time-consuming and costly, and the results are inherently dependent on the specific section of tissue that has been biopsied, which might not accurately reflect the entire tumor landscape and hinder, at least partially, the detection of MRD. Furthermore, new clones carrying additional variants can emerge in response to selective pressures, which may be potentially exacerbated by the administration of adjuvant therapy, entailing an increase of false negatives.

An interesting —although hardly explored— approach to partially circumvent these limitations would be to couple tumor-informed assays with a tumor-agnostic panel to account for uncharacterized or appearing mutations. One study assessed the benefit of adding a small tumor-agnostic panel of 10 CRC-related genes to a tumor-informed assay. In comparison to a conventional tumor-informed method, this combines approached raised the sensitivity of detection (96.0% vs. 88.0%) in a cohort of 34 stage I-IV CRC patients. More research is needed to validate this approach [@chenTumorinformedPatientspecificMRD2023].

[@chanTumorinformedTumoragnosticCirculating2023]

[Patient-specific tumor-informed circulating tumor DNA (ctDNA) analysis for molecular residual disease (MRD) detection in surgical patients with stage I-IV colorectal cancer (CRC). | Journal of Clinical Oncology](https://ascopubs.org/doi/abs/10.1200/JCO.2023.41.4_suppl.213)

### Discussion of studies comparing both approaches

Few studies have directly compared both approaches.

![[Tumor-informed versus tumor-agnostic approaches in addressing minimal residual disease — unraveling the optimal strategy (planning)#Studies comparing both approaches]]

## V. Emerging Technologies and Future Directions

- Fragmentomics
- CTCs
- Imaging techniques (?) (guided by AI?)

> 1. NGS -> Signatera and WES approach. brPROPHET
> 2. Epigenomics and beyond -> LUNAR-1, Guardant Health (ECLIPSE trial, failure of COBRA, etc.). Incorporation of other molecular markers cfRNA, miRNA...
> 3. Liquid biopsy -> CTCs
> 4. Imaging techniques, guided by AI.

[brPROPHET May Better Detect ctDNA and Identify MRD in Colorectal Cancer](https://www.targetedonc.com/view/brprophet-may-better-detect-ctdna-and-identify-mrd-in-colorectal-cancer)

### 1. Genomics

The advent of efficiency improvements in sequencing technologies and the development of new bioinformatic methods have enabled less costly assessments of greater genome regions.

#### Signatera assay: GALAXY trial and BESPOKE trial

The Signatera assay (Natera) has become a prominent advocate of the tumor-informed approach. In brief, initial tumor variants are sought after in both primary tumor tissue and matched normal tissue using a whole exome sequencing approach, which is a step beyond the targeted sequencing used by previous studies. A set of 16 or more mutations is worked out and subsequently considered for the creation of a bespoke test based on ultradeep multiplex PCR, which is then used to detect ctDNA and monitorize MRD [@SignateraOverview].

One key study that employs this approach in CRC is the GALAXY trial (UMIN000039205), part of the CIRCULATE-Japan study. Here, a cohort of 1039 patients with resectable CRC (stages II-IV) was subject to an observational study, with postoperative ctDNA status determination and follow-up to a median of 16.74 months. After surgery, the ctDNA status was used as the criterion for patient assignment to one of two related clinical trials: patients that tested positive for ctDNA were assigned to the ALTAIR trial, assessing the adequacy of chemotherapy escalation, and ctDNA-negative patients were assigned to the VEGA trial, studying the noninferiority of treatment deescalation. The ctDNA determination exhibited prognostic value, with patients that tested positive at 4 weeks after surgery showing an increased risk of recurrence (HR = 10.0) [@kotaniMolecularResidualDisease2023].

Signatera is currently also the centerpiece of the ongoing BESPOKE clinical trial (NCT04264702). This multi-center study tracks the ctDNA status of a prospective cohort of 2000 stage I-IV CRC patients, while assessing how this testing can affect the adjustment of adjuvant chemotherapy and patient outcomes. By using Signatera to detect the ctDNA, the results of this study may further contribute to its validation, and also to expand upon the current evidence of the clinical utility of the tumor-informed approach [@kasiBESPOKEStudyProtocol2021]. Recently presented preliminar data showed promising results in favor of the prognostic value of the test: in a subset of 350 patients (CRC stage II-III) of the BESPOKE cohort, patients with detectable minimal residual disease exhibited a reduced DFS (HR = 20.8; 95%CI: 10.0-43.4; p < 0.0001). Interestingly, while adjuvant chemotherapy proved to increase the median DFS of MRD-positive patients (18.7 vs 6.7 months; HR = 3.9; 95%CI: 1.3-11.5: p = 0.01), the same was not true for MRD-negative patients.

### 2. Epigenomics and beyond

> The advent of multi-omic integration most benefits tumor-agnostic assays, as they have shown in recent assays.

> As evidenced by the LUNAR-1 panel, tumor-agnostic assays gain the most benefits of the advent 


The study of whole genome methylation markers through techniques such as cfMeDIP-seq —which enables high resolution profiling at a relatively low cost when compared with bisulfite methods— is also starting to render insightful results. 

The study of cfDNA fragmentation patterns —a field termed **fragmentomics**— is also starting to gain traction in the setting of tumor-agnostic MRD detection. 

A study in 51 metastatic cancer patients from the INSPIRE trial (NCT02644369) that followed a tumor-agnostic approach showed that cancer specific methylation predicted OS and PFS similarly to mutation concentration in ctDNA. When integrating fragmentomics, both OS and PFS prediction surpassed that of mutation concentration [@zhao1664MOTumornaiveMethylomes2022].

The integration of multiple omics is certainly a target to pursue for the promise of a more sensitive and specific detection of ctDNA. [Redirecting](https://doi.org/10.1016/j.annonc.2022.07.1744). Similar approaches could be followed for MRD detection, although the high related costs would be currrently unaffordable for the daily clinical practice. 



EC:
- ECLIPSE
- GRAIL ??'
- MEDAL

The potential of 

### 2. Liquid Biopsy Technologies


Other liquid biopsy analytes beyond ctDNA have begun to be explored to address minimal residual disease, mostly in a tumor-agnostic fashion. A promising biomarker are circulating tumor cells (CTCs) —



- CTCs [@marcuelloCirculatingBiomarkersEarly2019]
- miRNA [@marcuelloCirculatingBiomarkersEarly2019]
- mRNA !!!!!!!!!!!!!!!!
- cfRNA (gtcS, [Fetching Title#tv5v](https://genomictestingcooperative.com/tumor-informed-vs-tumor-naive/))

The study of circulating tumor cells (CTCs) -cancer cells that are released from the tumor mass and migrate into preripheral circulation- is currently a hot topic of research.



### 4. Advanced Imaging Techniques

> The application of AI in radiomics.
> However, ctDNA still takes 


[@tibermacineRadiomicsModellingRectal2021]

### 4. AI

> Beyond radiomics, AI also has interesting applications...


