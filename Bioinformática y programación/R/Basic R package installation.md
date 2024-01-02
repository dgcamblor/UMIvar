---
tags:
  - how-to
---

Installation of basic R packages.
# Packages from CRAN

The package `languageserver` is a dependency for some tools, like the R extension in VSCode.

```
install.packages("languageserver")
```

```r
# Install basic packages
install.packages(c("tidyverse",
                   "devtools",
                   "rmarkdown"))

# Packages from Bioconductor

```r
# Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```