\name{Permanent-functions}
\alias{Permanent-functions}
\alias{cxPerm}
\alias{rePerm}
\alias{cxPermMinors}
\title{Functions for evaluating matrix permanents}
\description{
  These three functions are used in the classical Boson Sampling problem
}
\usage{
cxPerm(A)
rePerm(B)
cxPermMinors(C)
}
\arguments{
  \item{A}{a square complex matrix.}
  \item{B}{a square real matrix.}
  \item{C}{a rectangular complex matrix where \code{nrow(C) = ncol(C) + 1}.}
}
\value{
  \code{cxPerm(A)} returns a complex number: the permanent of the complex matrix \code{A}.\cr
  \code{rePerm(B)} returns a real number: the permanent of the real matrix \code{B}.\cr
  \code{cxPermMinors(C)} returns a complex vector of length \code{ncol(C)+1}: the permanents of all 
  \code{ncol(C)}-dimensional square matrices constructed by removing individual rows from \code{C}.
}
\details{
Permanents are evaluated using Glynn's formula (equivalently that of Nijenhuis and Wilf (1978))
}
\references{
  Glynn, D.G. (2010) The permanent of a square matrix. \emph{European Journal of Combinatorics}, \bold{31}(7):1887--1891.\cr
  Nijenhuis, A. and Wilf, H. S. (1978). \emph{Combinatorial algorithms: for computers and calculators.} Academic press.
}
\examples{
  set.seed(7)
  n <- 20
  A <- randomUnitary(n)
  cxPerm(A)
  #
  B <- Re(A)
  rePerm(B)
  #
  C <- A[,-n]
  v <- cxPermMinors(C)
  #
  # Check Laplace expansion by sub-permanents
  c(cxPerm(A),sum(v*A[,n]))
}

