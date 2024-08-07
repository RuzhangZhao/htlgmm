% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/htlgmm.R
\name{prod_fast}
\alias{prod_fast}
\title{Fast Matrix Dot Product with Transpose to Replace "\%*\%".}
\usage{
prod_fast(x, y)
}
\arguments{
\item{x}{The input matrix.}

\item{y}{The input matrix or vector or NULL.}
}
\value{
The dot product of 'x' and 'y'.
}
\description{
Fast Matrix Dot Product with Transpose to Replace "\%*\%".
}
\details{
Implemented by RcppEigen
}
