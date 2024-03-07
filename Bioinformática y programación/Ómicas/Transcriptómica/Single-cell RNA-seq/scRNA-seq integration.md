
The integration of scRNA-seq data can refer to integrating:

- Data from different experimental batches.
- From different donors.
- From different conditions.

A key example of scRNA-seq data integration is [@nietoSinglecellTumorImmune]. Seurat has an standard pipeline: [Introduction to scRNA-seq integration • Seurat](https://satijalab.org/seurat/articles/integration_introduction.html). This method is based on the identification of **anchor cells** between pairs of datasets, which are used to harmonize the integrated datasets. 

