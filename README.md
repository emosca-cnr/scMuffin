<img src="vignettes/images/scMuffin_logo_pontillism.jpg" width="200">

## scMuffin - A MUlti-Features INtegrative approach for SC data analysis

Single cell (SC) analysis is crucial to study the complex cellular heterogeneity of solid tumors, which is one of the main obstacles for the development of effective cancer treatments. Such tumors typically contain a mixture of cells with aberrant genomic and expression profiles affecting specific sub-populations that have a pivotal role in cancer progression, whose identification eludes bulk approaches. We present a MUlti-Features INtegrative approach for SC data analysis (scMuffin) that characterizes cell identity on the basis of multiple and complementary criteria. scMuffin provides functions to calculate a series of qualitative and quantitative scores, such as: expression of markers for normal and tumor conditions, pathway activity, cell hierarchy, multipotency state, copy number variations and cell cycle state. Cell-level scores are used for cell cluster annotation and combined to obtain alternative cell clusters. scMuffin integrates any type of cell- or cluster-associated data, and can be used for single-cell multi-omics analyses (e.g. mutations, gene expression). As a proof-of-principle, we studied a public dataset of human gliomas. scMuffin combines several tools to shed light on the identity of tumors cells and spot subtle cell types.

## Installation

scMuffin requires R >= 4.0.0 due to some of its dependencies, like Seurat [@Hao2021]. R can be installed from CRAN at the URL https://cran.r-project.org/index.html.

To succesfully install scMuffin you need some packages from Bioconductor (https://bioconductor.org) and github (https://github.com/). These packages can be installed using the following commands:

```{r, include=TRUE, eval=FALSE}
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(c("BiocStyle", "ComplexHeatmap", "DESeq2", "org.Hs.eg.db"))

if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("theislab/destiny")
```

The other dependencies, if missing, should be automatically installed using the following command:

```{r, include=TRUE, eval=FALSE}
devtools::install_github("emosca-cnr/scMuffin", build_vignettes = TRUE)
```

## Documentation
Please look at the vignette included in the R package:
```{r, eval=FALSE}
vignette("scMuffin")
```

