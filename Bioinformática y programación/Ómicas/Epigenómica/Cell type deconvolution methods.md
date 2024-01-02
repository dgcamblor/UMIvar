
## Epigenomic deconvolution

- CelFEER. [Cell type deconvolution of methylated cell-free DNA at the resolution of individual reads - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10236360/).
- **CancerLocator**. Determines the presence and location of tumors. Infers the proportion and tissue of origin of ctDNA in WGBS data [@kangCancerLocatorNoninvasiveCancer2017].
- **CancerDetector**. Identifies traces amounts of ctDNA in plasma (cfDNA), **at the level of individual reads** [@liCancerDetectorUltrasensitiveNoninvasive2018].  Works with sequencing data at **low to medium coverage (1x to 10x)**.

## Transcriptomic deconvolution

The main transcriptomic deconvolution methodologies are: 

-   [**Quantiseq**](https://icbi.i-med.ac.at/software/quantiseq/doc/). Deconvolución en 10 tipos inmunológicos, especialmente orientada a muestras tumorales. Trabaja en muestras de RNAseq. Puede identificar tipos celulares no caracterizados.
-   [**TIMER**](http://timer.cistrome.org/). Deconvolución de infiltrados inmunológicos en diferentes tipos de cáncer. En realidad, esta herramienta permite obtener abundancias estimadas por múltiples tipos de métodos de deconvolución, entre los que se encuentra
-   [**CIBERSORT/CIBERSORTx**](https://cibersortx.stanford.edu/). Diferenciación entre 22 células inmunes.
-   [**MCP-counter**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5). Método para la cuantificación de la abundancia absoluta de poblaciones del microambiente tumoral: 8 tipos inmunes y 2 estromales.
-   [**xCell**](https://xcell.ucsf.edu/). Contempla 64 tipos celulares inmunes y estromales.
-   [**EPIC**](http://epic.gfellerlab.org/). Deconvolución de tipos celulares inmunes y otras células no malignas en los tumores. Como punto extra, puede identificar tipos celulares no caracterizados.
-   [**ABIS**](https://giannimonaco.shinyapps.io/ABIS/). Orientado a muestras de células mononucleares en sangre periférica. 
-   [**ConsensusTME**](https://github.com/cansysbio/ConsensusTME). Deconvolución de células del microambiente tumoral.

A fundamental R package is [[immunedeconv]], which gathers all of these methodologies.