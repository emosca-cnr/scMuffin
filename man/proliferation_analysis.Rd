% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proliferation_analysis.R
\name{proliferation_analysis}
\alias{proliferation_analysis}
\title{Define a proliferation score}
\usage{
proliferation_analysis(
  scMuffinList = NULL,
  mc.cores = 1,
  nbins = 25,
  nmark_min = 5,
  ncells_min = 5,
  k = 99,
  kmin = 49,
  score_type = c("relative", "mean"),
  mean_scale = TRUE,
  gsl = NULL
)
}
\arguments{
\item{scMuffinList}{scMuffinList object}

\item{mc.cores}{number of cores}

\item{nbins}{number of bins to split the distribution of average gene expression}

\item{nmark_min}{number of minimum markers that are required for the succesful calculation of a signature}

\item{ncells_min}{number of minimum cells in which a gene set has to be succesfully calculated}

\item{k}{number of permutations}

\item{kmin}{minimum number of permutations; due to missing values it is hard to ensure that a signature can be compared to k permutations in every cell}

\item{score_type}{type of score. if "relative", then the score is the difference between the observed gene set average expression and that of a k permutations; if "mean" the score is equal to the observed gene set average expression}

\item{mean_scale}{whether to scale the values obtained using score_type="mean"}

\item{gsl}{list with two gene sets with 1/S and G2/M markers. If NULL gsls_Symbol$Tirosh is used. See gsls_Symbol$Tirosh to properly format custom gene sets.}
}
\value{
scMuffinList with element "proliferation", a list with:
\itemize{
\item{summary: data.frame with proliferation value}
\item{full: data.frame with the vales of each of the two considered gene sets}
}
}
\description{
Define a proliferation score on the basis of G1/S and G2/M markers (by default gsls_Symbol$Tirosh)
}
