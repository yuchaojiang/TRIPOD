% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/influence.R
\name{getCooksD}
\alias{getCooksD}
\title{Get Cook's distance}
\usage{
getCooksD(Yg, Xt, Yj, metacell.celltype, ...)
}
\arguments{
\item{Yg}{a numeric vector containing RNA expression values.}

\item{Xt}{a numeric vector containing chromatin accessibility.}

\item{Yj}{a numeric vector containing TF expression.}

\item{metacell.celltype}{a character vector specifying cell types of the metacells.}

\item{...}{further arguments to be passed.}
}
\description{
This function computes Cook's distance on the linear regression model
and prepares data for Cook's bar plot.
}
