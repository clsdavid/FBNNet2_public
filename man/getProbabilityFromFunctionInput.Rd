% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelling_FBN.R
\name{getProbabilityFromFunctionInput}
\alias{getProbabilityFromFunctionInput}
\title{An internal method to get the probability of regulatory function from the cube}
\usage{
getProbabilityFromFunctionInput(
  funcType,
  FBNExpression,
  FBNProbability,
  preGeneInputs
)
}
\arguments{
\item{funcType}{type of function}

\item{FBNExpression}{expression of function}

\item{FBNProbability}{The probability of a FBN connection}

\item{preGeneInputs}{pre input genes' state}
}
\value{
A probability of the target regulatory function
}
\description{
An internal method to get the probability of regulatory function from the cube
}
\examples{
##coming later
}
