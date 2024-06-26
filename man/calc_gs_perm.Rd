% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_gs_perm.R
\name{calc_gs_perm}
\alias{calc_gs_perm}
\title{Calculate permutations}
\usage{
calc_gs_perm(rll = NULL, perm = NULL, gs = NULL, fract_min = 0.2)
}
\arguments{
\item{rll}{numeric matrix of genes-by-ranking criteria; each column contains numeric values; rownames are mandatory}

\item{perm}{vector of permuted names}

\item{gs}{gene set}

\item{fract_min}{only cluster of size less or equal to this fraction of cell with not null feature values will be analysed}
}
\description{
Calculate permutations
}
