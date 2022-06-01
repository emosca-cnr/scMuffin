<img src="vignettes/images/scMuffin_logo_pontillism.jpg" width="200">

## scMuffin - A MUlti-Features INtegrative approach for SC data analysis

Single-cell (SC) gene expression analysis is crucial to dissect the complex cellular heterogeneity of solid tumors, which is one of the main obstacles for the development of effective cancer treatments. Such tumors typically contain a mixture of cells with aberrant genomic and transcriptomic profiles affecting specific sub-populations that might have a pivotal role in cancer progression, whose identification eludes bulk RNA-sequencing approaches. We present scMuffin, an R package that enables the characterization of cell identity in solid tumors on the basis of multiple and complementary criteria applied on SC gene expression data. scMuffin provides a series of functions to calculate several different qualitative and quantitative scores, such as: expression of marker sets for normal and tumor conditions, pathway activity, cell state trajectories, CNVs, chromatin state and proliferation state. Thus, scMuffin facilitates the combination of various evidences that can be used to distinguish normal and tumoral cells, define cell identities, cluster cells in different ways, link genomic aberrations to phenotypes and identify subtle differences between cell subtypes or cell states.

## Installation

scMuffin requires R >= 4.0.0 due to some of its dependencies, like Seurat (https://satijalab.org/seurat). R can be installed from CRAN at the URL https://cran.r-project.org/index.html.

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

