% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{getTFGenePairsForSingleGene}
\alias{getTFGenePairsForSingleGene}
\title{Get TF-gene pairs with significantly nonzero coefficients for a given target gene}
\usage{
getTFGenePairsForSingleGene(xymats, coef.name, pval.name, sign, fdr.thresh)
}
\arguments{
\item{xymats}{a list returned by {\code{\link{fitModel}}}.}

\item{coef.name}{a character string specifying the name of the coefficient.
This must be one of the names of the elements in the xymats object.}

\item{pval.name}{a character string specifying the name of the element in
the xymats object that stores the P value.}

\item{sign}{a character string specifying the sign.
This should be either "positive" or "negative".}

\item{fdr.thresh}{a numerical scalar setting a FDR threshold.}
}
\value{
A data frame with 8 columns containing the following:
\item{gene}{a character string containing the name of the target gene.}
\item{TF_num}{an integer representing a TF index.
The mapping between the indices and the TF names are the same as that
in the xymats object.}
\item{TF}{a chracter string representing the name of the peak.}
\item{coef}{a numerical scalar representing the value of the coefficient}
\item{pval}{a numerical scalar representing the nominal p value}
\item{adj}{a numerical scalar representing the BH-adjusted p value}
}
\description{
This function takes a list returned by {\code{\link{fitModel}}}
as input and generate a data frame containing TF-gene pairs with significantly
positive or negative coefficients for a given target gene.
}
