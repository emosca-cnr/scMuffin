% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_umap.R
\name{plot_umap}
\alias{plot_umap}
\title{plot_umap}
\usage{
plot_umap(
  Seu_obj,
  file = "umap.jpg",
  labels = NULL,
  group.by = NULL,
  feature_plot = FALSE,
  lab_size = 1,
  lab_color = "black",
  adj_outliers = FALSE,
  width = 180,
  height = 180,
  units = "mm",
  res = 300,
  text.size = 5,
  ...
)
}
\arguments{
\item{Seu_obj}{seurat object, object with saved dimension reduction components}

\item{file}{string, file name output}

\item{labels}{cluster labels}

\item{group.by}{a feature to colour by (e.g. cluster ID)}

\item{feature_plot}{whether to call Seurat::FeaturePlot}

\item{lab_size}{label size}

\item{lab_color}{label color}

\item{adj_outliers}{logical, whether to adjust the group.by scores, removing outliers}

\item{width}{image width}

\item{height}{image height}

\item{units}{image units}

\item{res}{image resolution}

\item{text.size}{text size}

\item{...}{further arguments for Seurat::FeaturePlot or Seurat::DimPlot}
}
\description{
Generate a UMAP visualization
}
\details{
Generate a UMAP visualization
}
