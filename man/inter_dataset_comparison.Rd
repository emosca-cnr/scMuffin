% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inter_datasets_comparison.R
\name{inter_dataset_comparison}
\alias{inter_dataset_comparison}
\title{Inter-dataset cluster similarity}
\usage{
inter_dataset_comparison(
  seu_obj_list = NULL,
  gsl = NULL,
  genes_min = 3,
  genes_max = 500,
  mc.cores = 1,
  ncells_min = 5,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  ...
)
}
\arguments{
\item{seu_obj_list}{list of Seurat objects}

\item{gsl}{gene set list}

\item{genes_min}{minimum number of genes required among the markers of a cluster}

\item{genes_max}{maximum number of genes required among the markers of a cluster}

\item{mc.cores}{number of cores}

\item{ncells_min}{minim number of cells in a cluster}

\item{cluster_rows}{whether to cluster or not the rows}

\item{cluster_columns}{whether to cluster or not the columns}

\item{...}{arguments passed to calculate_gs_scores}
}
\value{
A list with:
\itemize{
  \item{clust_sim, matrix of clucster similarity;}
  \item{m1=m1}
  \item{m2=m2}
  \item{markers_1, markers of dataset 1;}
  \item{markers_2, markers of dataset 2;}
}
}
\description{
Quantify the similarity between clusters of two datasets, on the basis of the average cluster marker expression
}
