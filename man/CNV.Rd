% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNV.R
\name{CNV}
\alias{CNV}
\title{CNV}
\usage{
CNV(x, wnd_size = 100, na.rm = FALSE)
}
\arguments{
\item{x}{list of dataframes retrieved by 'preprocess_object_for_cnv'.}

\item{wnd_size}{number of adjacent genes to consider}

\item{na.rm}{whether to remove genes that have missing value of expression or not.}
}
\description{
Function to be used in calculate_CNV.
}
\references{
"Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma" Patel et al. (DOI: 10.1126/science.1254257), Supplementary materials
}
