\name{bosonSampler}
\alias{bosonSampler}
\title{Function for independently sampling from the Boson Sampling distribution}
\description{
  The function implements the Boson Sampling algorithm defined in Clifford and Clifford (2017) \url{https://arxiv.org/abs/1706.01260}
}
\usage{
bosonSampler(A, sampleSize, perm = FALSE)
}
\arguments{
  \item{A}{the first \code{n} columns of an \code{(m x m)} random unitary matrix, 
  \cr see \code{\link{randomUnitary}}}
  \item{sampleSize}{the number of independent sample values required for given \code{A}}
  \item{perm}{TRUE if the permanents and pmfs of each sample value are required}
}
\details{
  Let the matrix \code{A} be the first \code{n} columns of an \code{(m x m)} random unitary matrix, then 
  \code{X <- bosonSampler(A, sampleSize = N, perm = TRUE)} provides \code{X$values}, \code{X$perms} and \code{X$pmfs}, 
  \cr
  \cr The component \code{X$values} is an \code{(n x N)} matrix with columns that are 
  independent sample values from the Boson Sampling distribution. 
  Each sample value is a vector of \code{n} integer-valued output modes in random order. The elements of the vector can be sorted in 
  increasing order to provide a multiset representation of the sample value. 
  \cr 
  \cr The outputs \code{X$perms} and \code{X$pmfs} are 
  vectors of the permanents and probability mass functions (pmfs) associated with the sample values.
  The permanent associated with a sample value \code{v = (v_1,...,v_n)} is the permanent of an \code{(n x n)} matrix constructed with rows 
  \code{v_1,...,v_n} of \code{A}. Note the constructed matrix, \code{M}, may have repeated rows since \code{v_1,...,v_n} 
  are not necessarily distinct.  
  The pmf is calculated as \code{Mod(pM)^2/prod(factorial(tabulate(c))} where \code{pM} is the permanent of \code{M}.
}
\value{
\code{X = bosonSampler(A, sampleSize = N, perm = TRUE)} provides \code{X$values}, \code{X$perms} and \code{X$pmfs}. See Details.
}
\references{
  Clifford, P. and Clifford, R. (2017) The Classical Complexity of Boson Sampling, \url{https://arxiv.org/abs/1706.01260}
}
\examples{
set.seed(7)
n <- 20  # number of photons
m <- 200 # number of output modes
A <- randomUnitary(m)[,1:n]
# sample of output vectors
valueList <- bosonSampler(A, sampleSize = 10)$values 
valueList
# sample of output multisets 
apply(valueList,1, sort) 
#
set.seed(7)
n <- 12  # number of photons
m <- 30 # number of output modes
A <- randomUnitary(m)[,1:n]
# sample of output vectors
valueList = bosonSampler(A, sampleSize = 1000)$values 
# Compare frequency of output modes at different
# positions in the output vectors 
matplot(1:m,apply(valueList,1,tabulate), pch =20, t = "p", 
xlab = "output modes", ylab = "frequency")
}