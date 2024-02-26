---
tags:
  - paper
  - phd
---
> Regardless of the technique, ctDNA assays for the detection of MRD will require high specificity to escalate adjuvant therapy in patients who wouldf be categorized as low risk according to clinicopathologic factors and would ultimately recur without additional therapy, and high sensitivity to allow safe de-escalation of therapy in those patients who are traditionally treated with intensive chemotherapy and avoidable adverse events. [@jacomeMinimalResidualDisease2023]. 
> Un poco a colación de lo que dicen Reinert 2024

>In a prospective study with 112 mCRC patients who had undergone liver resection with curative intent, ctDNA positivity by a tumor-informed assay was the most significant prognostic factor for disease-free survival (HR: 5.7, 95% CI 3.3–10.0) [52]. MRD was detected in 61 (54%) patients, of which 59 (97%) presented recurrence at the time of data cutoff.

## II. Tumor-agnostic strategies

> —also termed tumor-uninformed or tumor-naïve—

> In recently presented results, however, the [[COBRA trial]] was not able to meet its primary endpoint, and was promptly stopped due to futility in achieving statistical significance [@morrisPhaseIIResults2024].

 A prospective study on a cohort of 240 CRC patients stage II-III showed decrease 2-year RFS rates for ctDNA-positive patients (HR = 10.98, 95% CI = 5.41-22.72, p < 0.001), a difference that was further exacerbated when performing serial ctDNA analysis in a subcohort of 125 patients (24.0% vs. 96.0%; HR = 32.02, 95 CI = 10.79-95.08, p < 0.001). Importantly, this last analysis exhibited a sensitivity of 82.6% and a specificity of 94.1% for the recurrence detection [@chenPostoperativeCirculatingTumor2021].
### Methylation

[Efficient detection and post-surgical monitoring of colon cancer with a multi-marker DNA methylation liquid biopsy | PNAS](https://www.pnas.org/doi/full/10.1073/pnas.2017421118)

Taieb, J.; Taly, V.; Henriques, J.; Bourreau, C.; Mineur, L.; Bennouna, J.; Desrame, J.; Louvet, C.; Lepere, C.; Mabro, M.; et al. Prognostic value and relation with adjuvant treatment duration of ctDNA in stage III colon cancer: A post-hoc analysis of the PRODIGE-GERCOR IDEA-France trial. Clin. Cancer Res. 2021. [CrossRef] [PubMed]

The use of techniques that detect epigenomic signatures assumes that hypermethylation of promoter regions of some tumor suppressor genes is a common finding in CRC. Hypermethylation of WNT-inhibitor-factor-1 (WIF1) and neuropeptide Y (NPY) has been validated as a marker of both advanced and localized colon cancer [35,36].

[Response prediction by mutation- or methylation-specific detection of ctDNA dynamics in pretreated metastatic colorectal cancer - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10541738/)

## III. Tumor-informed strategies

![[Tumor-informed versus tumor-agnostic approaches in addressing minimal residual disease — unraveling the optimal strategy (planning)#Tumor-informed]]

| Study | Approach | Sensitivity | Specificity | Cohort |
| ---- | ---- | ---- | ---- | ---- |
|  |  |  |  |  |

- MEDOCC-CrEATE trial
	- Protocol: [Circulating tumor DNA guided adjuvant chemotherapy in stage II colon cancer (MEDOCC-CrEATE): study protocol for a trial within a cohort study | BMC Cancer | Full Text](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07252-y)
	- Currently enrolling: search for more

Throughout the literature, tumor-informed approaches have by far been the most frequently employed strategy to assess the presence of ctDNA an MRD.

### ddPCR 

One approximation extended among the initial studies that assessed ctDNA as a proxy for MRD detection was the use of droplet digital PCR (ddPCR) to track the presence of specific mutations previously identified in the tumor tissue. A pioneer study by Reinert et al. used a ddPCR-based approach to evaluate the presence of large structural variants found in the tumor tissue in a cohort of 118 patients. Although at the postoperative landmark analysis both the sensitivity and specificity reported for the prediction of recurrence were of 100%, they were based on 6 relapsing and 5-non relapsing patients, which made difficult to draw conclusions by itself. However, the analysis of CEA —the gold-standard biomarker for MRD detection— in the same patients rendered a sensitivity and a specificity of 67% and 100%, respectively, which hinted at a trend that would be later confirmed by larger studies [@reinertAnalysisCirculatingTumour2016].

Subsequent studies shifted towards focusing on the use of single nucleotide variants (SNVs) for ctDNA calling. In a prospective study of 150 patients with CRC localized disease, tumor tissue was assesed with a 29 gene custom panel, and the two most prevalent mutations were tracked using ddPCR. Post-operative ctDNA detection was significantly correlated with worse DFS as an independent variable (HR = 11.64, 95% CI = 3.67-36.88, p = 0.0001). Of note, while the landmark analysis detected ctDNA in 47.1% (of 17 patients) that finally experienced relapse, this proportion was increased to 82.4% when performing serial ctDNA analysis [@tarazonaTargetedNextgenerationSequencing2019].

The use of ddPCR for the detection of ctDNA is still a method of great relevance in current research due to its benefits in terms of sensitivity, cost and turnaround time. A recent prospective study in a large cohort of 851 CRC patients stage II-III performed whole-exome sequencing (WES) on tumor tissue and tracked patient specific mutations via ddPCR. Notably, ctDNA detection on the ddPCR data was based on the consensus of a set of four ctDNA-calling algorithms. At the landmark timepoint (a median of 21 days after surgery), this approach reported a sensitivity of 35% and a specificity of 98% for detecting reccurency [@henriksenUnravelingPotentialClinical2024]. Indeed, this is the method of choice for the ongoing prospective, randomized clinical trial IMPROVE-IT2, which aims at studying whether the ctDNA detection can guide the postoperative computed tomography surveillance to produce a clinical impact in the management of CRC stages II-III [@norsIMPROVEIT2ImplementingNoninvasive2020].

### SafeSeqS approach: the DYNAMIC studies

In contrast to ddPCR, the NGS-based assessment of ctDNA has captured most of the attention of the current literature. One of the earliest key studies that supported this approach was led by Tie et al., and was performed in a prospective cohort of 230 CRC stage II patients (trial ACTRN12612000326897), where postoperative ctDNA-positive status was related to higher risk of recurrence among those patients that were not treated with adjuvant chemotherapy (HR = 18; 95%CI = 7.9-40). To discern the presence of ctDNA/MRD, the group initially studied FFPE tumor tissues in search for the variant with the highest allele frequency, using a targeted multiplex PCR aiming for 15 genes relevant to CRC cancer [@tieCirculatingTumorDNA2016]. Subsequently, those predominant variants were tracked using personalized Safe-SeqS assays, which rely on the use of Unique Molecular Identifiers (UMIs) to enable enhanced error correction and calling of low frequency variants [@kindeDetectionQuantificationRare2011]. The finding of said mutations was used as the criterion for patient classification into ctDNA-positive or -negative groups. This initial approximation rendered a sensitivity of 48% and a specificity of 100% for predicting recurrence at a timepoint of 36 months, which pointed to the great capacity of tumor-informed assays to discard MRD and hinted at the need of tracking multiple variants in each patient for increasing sensitivity [@tieCirculatingTumorDNA2016]. The same approach was used to describe a similar trend in RFS reduction for ctDNA-positive patients in stage III CRC (n = 96; HR = 3.8, 95% CI: 2.4-21.0, p < 0.001) and locally advanced rectal cancer (LARC) (n = 159; HR = 13.0, 95% CI: 5.5-31, p < 0.001), which predicted recurrence in 42% (10/24) and 58% (11/19) patients, respectively [@tieCirculatingTumorDNA2019; @tieSerialCirculatingTumour2019].

The growing evidence for the clinical utility of ctDNA in the postoperative setting encouraged the same group to conduct a series of related clinical trials aimed at evidencing the benefit of a ctDNA-guided handling of the adjuvant therapy, comprised of treatment escalation in ctDNA-positive patients and de-escalation in ctDNA-negative patients. These prospective interventional studies —DYNAMIC, DYNAMIC-III and DYNAMIC-Rectal— assess patients with CRC stage II, stage III and locally advanced rectal cancer, respectively, which undergo randomization after surgery into one of two groups: patients receiving treatment according to their ctDNA status and patients receiving conventional treatment. Importantly, to determine the presence of ctDNA, the studies adopted the aforementioned Safe-SeqS tumor-informed methodology.

The DYNAMIC trial (ACTRN12615000381583) was the first to reveal a clinical benefit of the ctDNA-guided management of adjuvant therapy in colon cancer stage II. By random assignment, patients were assigned to an standard management group (n = 153) or a ctDNA-guided management group (n = 302). In this approximation, patients that were ctDNA-positive at 4 or 7 weeks after surgery were recommended to receive adjuvant oxaliplatin- or fluoropyrimidine-based chemotherapy. Conversely, patients testing negative for ctDNA did not receive adjuvant chemotherapy. In the published results, RFS for the ctDNA-guided group was comparable to that of the standard management group (93.5% vs. 92.4%), indicating the suitability of chemotherapy reduction in ctDNA negative stage II patients [@tieCirculatingTumorDNA2022]. The DYNAMIC-III and DYNAMIC-RECTAL trials have yet to report their results.

> [Program Guide – ASCO Meeting Program Guide](https://meetings.asco.org/abstracts-presentations/228859) 
> [Studies Confirm the Utility of ctDNA in Guiding Adjuvant Chemotherapy in CRC | ASCO Daily News](https://dailynews.ascopubs.org/do/studies-confirm-utility-ctdna-guiding-adjuvant-chemotherapy-crc)

### WES in tissue characterization: Signatera

> Each target can achieve an average depth of 100,000x with ultra‑deep sequencing, and this approach has exhibited a sensitivity >98% at 0.01‑0.02% ctDNA concentrations (60). [@chenDetectingLiquidRemnants2023] The predefined threshold is based on Natera’s proprietary variant calling method wherein detecting at least 2 out of 16 variants ensures the optimal analytical performance of the assay with > 95% sensitivity at 0.01% mean variant allele frequency and with 99.7% specificity [@kotaniMolecularResidualDisease2023].

The commercial assay Signatera (Natera) has currently become a prominent advocate of the tumor-informed approach, which relies on WES and multiplex PCR for detecting the main 16 tumor-specific mutations. A prospective study by Reinert et al., in a total of 130 CRC patients stage I-III reported a significantly high RFS for the positive status of ctDNA both at the timepoint of 30 days (HR = 7.2, 95% CI: 2. -19.0, p < .001) and, markedly, in the serial ctDNA analysis (HR = 43.5, 95% CI, 9.8-193.5 p < .001). Importantly, in this surveillance setting, the reported values  for sensitivity and specificity were 88% and 98%, respectively [@reinertAnalysisPlasmaCellFree2019]. Similarly, in another study on a cohort of 168 stage III patients, ctDNA determination using Signatera was able to predict relapse with a sensitivity of 42% (16/38), while this proportion increased to 88% (21/24) in the serial analysis [@henriksenCirculatingTumorDNA2022].

Based upon the previous studies that reported results with Signatera, large-scale clinical trials employing this approach in CRC have emerged, such as the GALAXY trial (UMIN000039205). This is an observational part of the CIRCULATE-Japan study, and it aims to assess the presence of ctDNA in the postoperative setting for a subcohort of stage II-IV patients. Importantly, the observed ctDNA status is used in the framework of the CIRCULATE-Japan as the criterion for patient assignment to one of two related clinical trials: patients that test positive for ctDNA are assigned to the ALTAIR trial, assessing the adequacy of chemotherapy escalation, and ctDNA-negative patients are assigned to the VEGA trial, studying the non-inferiority of treatment de-escalation [@taniguchiCIRCULATEJapanCirculatingTumor2021]. The first results reported from the GALAXY trial came from an interim analysis on a cohort of 1039 patients spanning the aforementioned stages II-IV. Patients were followed-up postoperatively for up to a median of 16.74 months and, at the landmark timepoint of 4 weeks after surgery, the relapse rate for ctDNA-negative patients was 9.5% (of 852), while for ctDNA-positive patients amounted to 61.4% (of 187) (HR = 10.0; 95% CI = 7.7 - 14.0, p < 0.0001) [@kotaniMolecularResidualDisease2023].

Signatera is currently also the centerpiece of the ongoing BESPOKE clinical trial (NCT04264702). This multi-center study aims to track the ctDNA status of a prospective cohort of 2000 stage I-IV CRC patients, while assessing how this testing can affect the adjustment of adjuvant chemotherapy and patient outcomes. By using Signatera to detect the ctDNA, the results of this study may further contribute to its validation, and also to expand upon the current evidence of the clinical utility of the tumor-informed approach [@kasiBESPOKEStudyProtocol2021]. Recently presented preliminar data showed promising results in favor of the prognostic value of the test: in a subset of 350 patients (CRC stage II-III) of the BESPOKE cohort, patients with detectable minimal residual disease exhibited a reduced DFS (HR = 20.8; 95%CI: 10.0-43.4; p < 0.0001). Interestingly, while adjuvant chemotherapy proved to increase the median DFS of MRD-positive patients (18.7 vs 6.7 months; HR = 3.9; 95%CI: 1.3-11.5: p = 0.01), the same was not true for MRD-negative patients [@kasiCirculatingTumorDNA2024].

> Another similar technology is AlphaLiquid Detect [Personalised circulating tumour DNA assay with large-scale mutation coverage for sensitive minimal residual disease detection in colorectal cancer | British Journal of Cancer](https://www.nature.com/articles/s41416-023-02300-3).

[FoundationOne Tracker | Foundation Medicine](https://www.foundationmedicine.com/test/foundationone-tracker) [Study Details | | ClinicalTrials.gov](https://clinicaltrials.gov/study/NCT04259944)

## IV. Comparison

### Strengths and weaknesses

As shown by the literature, tumor-informed approaches can detect MRD with higher sensitivity and specificity by looking for specific known mutations. In simple terms, this is possible because the search is narrowed down to allow for more coverage in informative regions, rather than relying on exploring a panel of genes with less coverage, and which may or may not be of enough relevance for the patient's tumor. This method also ensures that germline variants and variants related to clonal hematopoiesis of indeterminate potential are directly filtered out of the analysis [@steensmaClonalHematopoiesisIndeterminate2015]. Although the same can be achieved in tumor-agnostic methods by also performing paired sequencing in peripheral blood, in the case of CHIP mutations some of them may be unique to cfDNA, causing false positives that will not be sought after in tumor-informed approaches.

Tumor-informed assays are, however, inherently limited by their reliance on tumor tissue. The most obvious examples in which this approximation is unfeasible are the cases where the tumor tissue biopsy is not possible, either due to inaccessibility or to the available sample not being adequate for genomic characterization, such as after extensive neoadjuvant treatment. Whenever the tumor tissue is available, the tissue biopsy and subsequent sequencing for the creation of a custom assay are time-consuming and costly, potentially impacting their clinical adequacy. The results are also inherently dependent on the specific section of tissue that has been biopsied, which might not accurately reflect the entire tumor landscape and hinder, at least partially, the detection of MRD. Furthermore, new clones carrying additional variants can emerge in response to selective pressures —a process exacerbated by the administration of adjuvant therapy— entailing an increase of false negatives.

An interesting approach —although hardly explored in CRC— to partially circumvent these limitations would be to couple tumor-informed assays with a tumor-agnostic panel to account for uncharacterized or appearing mutations. One study assessed the benefit of adding a small tumor-agnostic panel of 10 CRC-related genes to a tumor-informed assay. In comparison to a conventional tumor-informed method, this combines approached raised the sensitivity of ctDNA detection (96.0% vs. 88.0%) in a cohort of 34 stage I-IV CRC patients [@chenTumorinformedPatientspecificMRD2023]. More research is needed to explore the benefits of this approach.

[@chanTumorinformedTumoragnosticCirculating2023]

### Discussion of studies comparing both approaches

Few studies have directly compared both approaches.


## V. Emerging Technologies and Future Directions

### 1. Breakthroughs in ctDNA genomics

With NGS having become the most prevalent method for the molecular characterization of tumor tissue and ctDNA, the advent of efficiency improvements in sequencing technologies and the development of new bioinformatic methods have revolutionized both tumor-informed and tumor-agnostic strategies by enabling less costly assessments of greater genome regions, and facilitating the implementation into clinical practice.

#### WES for better tumor tissue characterization

Whole exome sequencing (WES) is becoming a established method for variant characterization; this approach aims at capturing the entire coding regions of the genome, which are essentially the most relevant for the tumor landscape. Thus, WES currently represents a good trade-off between the breadth of the genomic regions analyzed and the cost of sequencing them at a suitable coverture for detecting low frequency variants. 

The Signatera (Natera) commercial assay, which is currently offered as a lab service, best exemplifies the harnessing of NGS developments on the tumor characterization end. In brief, initial tumor variants are sought after in both primary tumor tissue and matched normal tissue using a WES approach. This broad approach allows to work out a set of 16 or more mutations that are subsequently considered for the creation of a bespoke test based on ultradeep multiplex PCR, with average coverages reaching values superior to 100,000X and variants being detected even at allele frequencies bellow 0.01% [@coombesPersonalizedDetectionCirculating2019]. Notably, the results of this test are processed via a propietary software for variant calling, considering the concomitant presence of at least 2 variants for ctDNA calling [@SignateraOverview]. The current evidence of its clinical utility has earned it the designation of Breakthrough Device by the Food and Drug Administration in the United States. As described before, many large-scale prospective clinical trials are being currently conducted that use Signatera as their method of choice for MRD-ctDNA detection (GALAXY, BESPOKE CRC), which will provide further evidence of the clinical utility of the method.

The commercial assay brPROPHET has also leveraged WES for tumor tissue characterization and variant targeting in ctDNA. In a preliminar report by Cao et al., up to 55 variants were tracked in the ctDNA of a cohort of 117 patients spanning stages I-IV. Although the conclusions on RFS and capability of predicting recurrence could not be draw due to the limited observation period, a direct comparison of the method with a fixed panel of 168 genes in the preoperative setting was provided. While brPROPHET detected ctDNA in 97.3% of the patients, the tumor-agnostic approach did so in 68.9%, highlighting the increased sensitivity of the tumor-informed approach. Further results that show the capacity of this assay remain to be seen [@caoPatientspecificTumorinformedCirculating2023].


> NGS approaches have become prevalent for tumor sequencing and have also been applied to ctDNA detection. WGS applied to cfDNA achieves a sequencing depth of 0.1× and WES achieves a sequencing depth of 100× ([Farris and Trimarchi, 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B28); [Heitzer et al., 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B46); [Murtaza et al., 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B67); [Cohen et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B19)). Although some studies suggest that WGS is feasible for clinical application to certain patients, it is prohibitive for routine clinical implementation of WGS because of its cost and time required to perform WGS and the associated bioinformatic analysis ([Welch et al., 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B101); [Chan et al., 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B10)). Therefore, WES turns out to be feasible to improve detection sensitivity and reduce cost while maintaining comprehensive coverage of likely mutated genomic regions. The exons are enriched for most of the pathogenic somatic mutations while they represent only 1.5% of the whole genome ([Choi et al., 2009](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448395/#B17)). Above all, there is an inverse correlation between sequencing breadth and detection cost, and sequencing depth _versus_ detection limit of detection. Due to the low level of ctDNA in body fluids, targeted approaches, including hybridization capture-based NGS, PCR amplicon-based NGS, are superior to more extensive sequencing approaches such as whole exome or whole genome sequencing.

#### WES in ctDNA characterization

On the ctDNA end of tumor characterization, the use of WES is also beginning to be explored in the context of MRD detection. Studies point to a great concordance of WES in circulating tumor DNA an tumor tissue, while better attending to intratumoral heterogeneity [@leenanitikulConcordanceWholeExome2023].

> Cite more papers
#### WGS

Whole-genome sequencing (WGS) has also began to be explored in tumor-informed approaches. MRD-Detect is a... [@zviranGenomewideCellfreeDNA2020]

The implementation of WGS in the clinical setting is, however, currently hindered by its high associated costs and the computational burden of the bioinformatic analysis. 

[Un momento…](https://www.researchgate.net/publication/361320611_Abstract_1959_Sensitive_detection_of_circulating_tumor_DNA_by_whole_genome_sequencing_Validation_of_MRDetect_using_serial_blood_samples_from_stage_III_colorectal_cancer_patients) -> Claus Lindbjerg, MRDdetect

### 2. Epigenomics and beyond

> 2. Epigenomics and beyond -> LUNAR-1, Guardant Health (ECLIPSE trial, failure of COBRA, etc.). Incorporation of other molecular markers cfRNA, miRNA...

Although genetic variants are by far the most assessed marker for detecting ctDNA/MRD, an increasing part of the attention has shifted to exploring the epigenomic landscape of ctDNA and, in particular, methylation. Indeed, aberrant methylation is a major player in tumorigenesis that, in contrast to tumor mutations, appears to be more consistent across tumors from different patients and stages [@garrigouStudyHypermethylatedCirculating2016; @haoDNAMethylationMarkers2017; @luoCirculatingTumorDNA2020]. As it can be inferred, the study of the methylome is promising for the agnostic assessment of MRD, in which sensitivity of detection is often hindered by the heterogeneity of the tumor mutational landscape.

The fragmentation patterns of the cfDNA (fragmentomics) also reveal important epigenomic insights that give away the presence of ctDNA: mainly, cancer-derived fragments are shorter than those of non-cancerous origin, and the specific regions in which the DNA is cleaved can show alterations in nucleosome positioning and chromatin accesibility. A key study in the field of fragmentomics was performed by Cristiano et al., which developed a gradient boosting machine learning model that was able to predict the presence of ctDNA by training on low depth WGS data from 236 patients with different cancer types and 245 healthy controls. The model was able to score an overall AUC of 0.94, and in CRC patients in particular, the sensitivity was of 88% at a specificity of 95% [@cristianoGenomewideCellfreeDNA2019].

Following the compelling evidence of the utility of epigenomic markers in the detection of ctDNA, the commercial assay Guardant Reveal (Guardant Health) has been developed to integrate the study of methylation and fragmentomics with the detection of genetic variants. 

Although the specific regions have not been disclosed, 


> Talk about Guardant Reveal and Guardant Infinity, which presents results in ASCO 2024, COSMOS study. ### Guardant Health to present data at ASCO GI supporting use of liquid biopsy to predict colon cancer recurrence [Guardant Health, Inc. - Guardant Health to present data at ASCO GI supporting use of liquid biopsy to predict colon cancer recurrence](https://investors.guardanthealth.com/press-releases/press-releases/2024/Guardant-Health-to-present-data-at-ASCO-GI-supporting-use-of-liquid-biopsy-to-predict-colon-cancer-recurrence/default.aspx)

> The advent of multi-omic integration most benefits tumor-agnostic assays, as they have shown in recent assays.

> As evidenced by the LUNAR-1 panel, tumor-agnostic assays gain the most benefits of the advent 

> The integration of genomics and epigenomics is best exemplified with guardant reveal, guardant infinity.

> WGBS is starting to be explored in the context of MRD detection (?) Although the relatively high cost makes its implementation

The study of whole genome methylation markers through techniques such as cfMeDIP-seq —which enables high resolution profiling at a relatively low cost when compared with bisulfite methods— is also starting to render insightful results. 

A study in 51 metastatic cancer patients from the INSPIRE trial (NCT02644369) that followed a tumor-agnostic approach showed that cancer specific methylation predicted OS and PFS similarly to mutation concentration in ctDNA. When integrating fragmentomics, both OS and PFS prediction surpassed that of mutation concentration [@zhao1664MOTumornaiveMethylomes2022].

The integration of multiple omics is certainly a target to pursue for the promise of a more sensitive and specific detection of ctDNA. [Redirecting](https://doi.org/10.1016/j.annonc.2022.07.1744). Similar approaches could be followed for MRD detection, although the high related costs would be currrently unaffordable for the daily clinical practice. 

> Discutir el papel de montagut [High accuracy of a blood ctDNA-based multimodal test to detect colorectal cancer - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0923753423040103#bib34)

EC:
- ECLIPSE
- MEDAL

The potential of 

### 3. Potential new biomarkers in liquid biopsy

- CTCs [@marcuelloCirculatingBiomarkersEarly2019]
- miRNA [@marcuelloCirculatingBiomarkersEarly2019]
- mRNA !!!!!!!!!!!!!!!!
- cfRNA (gtcS, [Fetching Title#tv5v](https://genomictestingcooperative.com/tumor-informed-vs-tumor-naive/))

Other liquid biopsy analytes beyond ctDNA have begun to be explored to address minimal residual disease, mostly in a tumor-agnostic fashion. 

Circulating tumor cells are also a potential MRD biomarker that has been explored in the tumor-agnostic setting. These are cells that have shed from the tumor and migrate into the peripheral blood circulation. Research has shown that the detection of CTCs in the postoperative setting is correlated with worse RFS in CRC [@rahbariMetaanalysisShowsThat2010]. Indeed, a prospective study in 130 CRC patients stage II-III hinted at the clinical utility of their use as an independent marker for predicting recurrence (HR = 2.739, 95% CI=1.042 - 7.200, p = 0.041) [@wangPrognosticModelsBased2019]. However, the study of CTCs is currently hindered by a lack of standardized methods for their detection and characterization, leading to limitations in the sensitivity and specificity of this approach.




### 4. Advanced Imaging Techniques

> The application of AI in radiomics.
> However, ctDNA still takes 


[@tibermacineRadiomicsModellingRectal2021]

### 4. (Bioinformatic methods and) AI

> Beyond radiomics, AI also has interesting applications...
> (CASTLE,19 Poisson,19,20 ALPACA,19,21 and Dynamic LOB19) [@henriksenUnravelingPotentialClinical2024]

## Challenges and future directions

> False negatives and positives, landmark vs sequential...

> Importance of the algorithm of ctDNA calling. Which is the cutoff for saying a patient is ctDNA positive?
