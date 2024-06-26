% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transcr_compl.R
\name{transcr_compl}
\alias{transcr_compl}
\title{Transcriptional Complexity and Entropy}
\usage{
transcr_compl(
  scMuffinList = NULL,
  min_counts = 5,
  min_cells = 10,
  min_genes = 500
)
}
\arguments{
\item{scMuffinList}{scMuffinList object}

\item{min_counts}{minimum number of counts}

\item{min_cells}{minimum number of cells in which a gene must have at least min_counts counts}

\item{min_genes}{minimum number of genes that a cell must express with at least min_counts counts}
}
\value{
scMuffinList with transcr_compl element a list with summary (full is empty):
\itemize{
  \item{tot_counts, total number of transcripts;}
  \item{n_genes, total number of expressed genes;}
  \item{C, transcriptional complexity;}
  \item{H, transcriptional entropy;}
}
}
\description{
Transcriptional Complexity and Entropy
}
