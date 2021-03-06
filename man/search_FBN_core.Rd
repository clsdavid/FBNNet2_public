% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_FBN.R
\name{search_FBN_core}
\alias{search_FBN_core}
\title{Internal method}
\usage{
search_FBN_core(
  fbnGeneCube,
  genes,
  useParallel = FALSE,
  threshold_confidence = 1,
  threshold_error = 0,
  threshold_support = 1e-04,
  maxFBNRules = 5
)
}
\arguments{
\item{fbnGeneCube}{A pre constructed Orchard cube}

\item{genes}{The target genes in the output}

\item{useParallel}{An option turns on parallel}

\item{threshold_confidence}{A threshold of confidence (between 0 and 1)
that used to filter the Fundamental Boolean functions}

\item{threshold_error}{A threshold of error rate (between 0 and 1) that
used to filter the Fundamental Boolean functions}

\item{threshold_support}{A threshold of support (between 0 and 1) that
used to filter the Fundamental Boolean functions}

\item{maxFBNRules}{The maximum rules per type (Activation and Inhibition)
per gene can be mined or filtered, the rest will be discarded}
}
\description{
Internal method
}
