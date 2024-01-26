---
tags:
  - paper
  - phd
---
> In recently presented results, however, the [[COBRA trial]] was not able to meet its primary endpoint, and was promptly stopped due to futility in achieving statistical significance [@morrisPhaseIIResults2024].
## II. Tumor-agnostic strategies

### Methylation


#### Genome-wide

[Response prediction by mutation- or methylation-specific detection of ctDNA dynamics in pretreated metastatic colorectal cancer - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10541738/)

## III. Tumor-informed strategies

- [ ] MEDOCC-CrEATE trial
	- Protocol: [Circulating tumor DNA guided adjuvant chemotherapy in stage II colon cancer (MEDOCC-CrEATE): study protocol for a trial within a cohort study | BMC Cancer | Full Text](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07252-y)
	- Currently enrolling: search for more

### DYNAMIC trial (SafeSeqS approach)

> Añadir sobre la capacidad de detección.
> El primero de TIE van solo a por una y el segundo por al menos una??? revisar

A growing body of research and clinical trial results has revealed the usefulness of tumor-informed approaches to MRD detection. 

One of the earliest key studies that supported this approach was performed in a prospective cohort of 230 CRC stage II patients (trial ACTRN12612000326897), where ctDNA detection after resection showed higher risk of recurrence as determined by a CT scan. Initially, the variant with the highest allele frequency in FFPE tumor tissue was determined for each patient using targeted multiplex PCR aiming for 15 genes relevant to CRC cancer. Afterwards, those variants were tracked using a personalized Safe-SeqS assay, an UMI-based approach that enables deep sequencing with error correction [@kindeDetectionQuantificationRare2011]. Patients were then classified as ctDNA-positive or -negative, depending on the finding of said mutations in plasma. The results showed that postoperative ctDNA-positive status was related to higher risk of recurrence among patients not treated with adjuvant chemotherapy (HR = 18; 95%CI = 7.9-40) [@tieCirculatingTumorDNA2016].

The same tumor-informed approach was adopted in the DYNAMIC trial (ACTRN12615000381583), which further studied the possibility of using the postoperative ctDNA status as a basis for making treatment decisions for patients with CRC stage II. By random assignment, patients were assigned to an standard management group (n = 153) or a ctDNA-guided management group (n = 302). In this approximation, patients that were ctDNA-positive at 4 or 7 weeks after surgery were recommended to receive adjuvant oxaliplatin- or fluoropyrimidine-based chemotherapy. Conversely, patients testing negative for ctDNA did not receive adjuvant chemotherapy. The RFS for the ctDNA-guided group was comparable to that of the standard management group (93.5% vs. 92.4%), indicating the suitability of chemotherapy reduction in ctDNA negative patients [@tieCirculatingTumorDNA2022].

> DYNAMIC-II es el que se ha comentado anteriormente. Está también el DYNAMIC-III (ACTRN12617001566325) en pacientes de CCR III, y el DYNAMIC-RECTAL (LARC?).

### Other studies

Studies using non-commercial tumor-informed assays have also rendered insightful results.

- [@wangPrognosticPotentialCirculating2019]
- One study used a custom targeted panel of 29 genes to study tumor tissue mutations, and performed a ctDNA-based follow-up of said variants using ddPCR. The detection of ctDNA after surgery was significantly associated to lower disease free survival (HR = 6.96), and so was the tracking across multiple postsurgical samples (HR = 8.03) [@tarazonaTargetedNextgenerationSequencing2019].
- [@reinertAnalysisPlasmaCellFree2019]
- [@chenPostoperativeCirculatingTumor2021]
- [@henriksenCirculatingTumorDNA2022]

Other studies that have since relied on custom

### Signatera assay: GALAXY trial and BESPOKE trial

The Signatera assay (Natera) has become a prominent advocate of the tumor-informed approach. In brief, initial tumor variants are sought after in both primary tumor tissue and matched normal tissue using a whole exome sequencing approach, which is a step beyond the targeted sequencing used by previous studies. A set of 16 or more mutations is worked out and subsequently considered for the creation of a bespoke test based on ultradeep multiplex PCR, which is then used to detect ctDNA and monitorize MRD [@SignateraOverview].

One key study that employs this approach in CRC is the GALAXY trial (UMIN000039205), part of the CIRCULATE-Japan study. Here, a cohort of 1039 patients with resectable CRC (stages II-IV) was subject to an observational study, with postoperative ctDNA status determination and follow-up to a median of 16.74 months. After surgery, the ctDNA status was used as the criterion for patient assignment to one of two related clinical trials: patients that tested positive for ctDNA were assigned to the ALTAIR trial, assessing the adequacy of chemotherapy escalation, and ctDNA-negative patients were assigned to the VEGA trial, studying the noninferiority of treatment deescalation. The ctDNA determination exhibited prognostic value, with patients that tested positive at 4 weeks after surgery showing an increased risk of recurrence (HR = 10.0) [@kotaniMolecularResidualDisease2023].

Signatera is also the centerpiece of the ongoing BESPOKE clinical trial (NCT04264702). This multi-center study tracks the ctDNA status of a prospective cohort of 2000 stage I-IV CRC patients, while assessing how this testing can affect the adjustment of adjuvant chemotherapy and patient outcomes. By using Signatera to detect the ctDNA, the results of this study may further contribute to its validation, and also to expand upon the current evidence of the clinical utility of the tumor-informed approach [@kasiBESPOKEStudyProtocol2021]. Recently presented preliminar data showed promising results in favor of the prognostic value of the test: in a subset of 350 patients (CRC stage II-III) of the BESPOKE cohort, patients with detectable minimal residual disease exhibited a reduced DFS (HR = 20.8; 95%CI: 10.0-43.4; p < 0.0001). Interestingly, while adjuvant chemotherapy proved to increase the median DFS of MRD-positive patients (18.7 vs 6.7 months; HR = 3.9; 95%CI: 1.3-11.5: p = 0.01), the same was not true for MRD-negative patients.

[Studies Confirm the Utility of ctDNA in Guiding Adjuvant Chemotherapy in CRC | ASCO Daily News](https://dailynews.ascopubs.org/do/studies-confirm-utility-ctdna-guiding-adjuvant-chemotherapy-crc)

## IV. Comparison

### Strengths and weaknesses

- Dependence on tumor tissue.
	- Problem if there is no tissue available. For example, in patients that have undergone neoadjuvant therapy, the resected specimens may have insufficient tissue or tumor content for genomic profiling [@chanTumorinformedTumoragnosticCirculating2023].
	- The requirement of tumor tissue comes at a evident cost, both in terms of time and economic resources.
	- The dependence on the initial characterization of the tumor tissue also raises the critical issue of missing on uncharacterized or appearing mutations. The mutations that are to be tracked are inherently dependent on the specific section of tissue that has been biopsied, which might not accurately reflect the entire tumor landscape. Furthermore, new clones carrying additional variants can emerge in response to selective pressures, which may be potentially exacerbated by the administration of adjuvant therapy. Both paths would entail an increase of false negatives in tumor-agnostic assays.
- 


[CHIP]



By sequencing both the tumor tissue and the ctDNA, variants related to clonal hematopoiesis of indeterminate potential can be directly filtered [@steensmaClonalHematopoiesisIndeterminate2015]. 

### Discussion of studies comparing both approaches

Few studies have directly compared both approaches.

![[Tumor-informed versus tumor-agnostic approaches in addressing minimal residual disease — unraveling the optimal strategy (planning)#Studies comparing both approaches]]

## V. Emerging Technologies and Future Directions

- Fragmentomics
- CTCs
- Imaging techniques (?) (guided by AI?)

### 1. Epigenomics and beyond

> The advent of multi-omic integration most benefits tumor-agnostic assays, as they have shown in recent assays.




The study of whole genome methylation markers through techniques such as cfMeDIP-seq —which enables high resolution profiling at a relatively low cost when compared with bisulfite methods— is also starting to render insightful results. 

The study of cfDNA fragmentation patterns —a field termed **fragmentomics**— is also starting to gain traction in the setting of tumor-agnostic MRD detection. 

A study in 51 metastatic cancer patients from the INSPIRE trial (NCT02644369) that followed a tumor-agnostic approach showed that cancer specific methylation predicted OS and PFS similarly to mutation concentration in ctDNA. When integrating fragmentomics, both OS and PFS prediction surpassed that of mutation concentration [@zhao1664MOTumornaiveMethylomes2022].

The integration of multiple omics is certainly a target to pursue for the promise of a more sensitive and specific detection of ctDNA. [Redirecting](https://doi.org/10.1016/j.annonc.2022.07.1744). Similar approaches could be followed for MRD detection, although the high related costs would be currrently unaffordable for the daily clinical practice. 


EC:
- ECLIPSE
- GRAIL ??'
- MEDAL

The potential of 

### 2. Liquid Biopsy Technologies: CTCs

- CTCs [@marcuelloCirculatingBiomarkersEarly2019]
- miRNA [@marcuelloCirculatingBiomarkersEarly2019]
- mRNA !!!!!!!!!!!!!!!!

The study of circulating tumor cells (CTCs) -cancer cells that are released from the tumor mass and migrate into preripheral circulation- is currently a hot topic of research.



### 3. Advanced Imaging Techniques

> The application of AI in radiomics.
> However, ctDNA still takes 


[@tibermacineRadiomicsModellingRectal2021]

### 4. AI

> Beyond radiomics, AI also has interesting applications...


