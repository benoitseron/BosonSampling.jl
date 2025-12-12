\name{randomUnitary}
\alias{randomUnitary}
\title{Random unitary}
\description{
Returns a square complex matrix sampled from the Haar random unitary distribution.
}
\usage{
randomUnitary(size)
}
\arguments{
\item{size}{dimension of matrix}
}
\value{
  A square complex matrix.
}
\examples{
m <- 25 # size of matrix (m x m)
set.seed(7)
U <- randomUnitary(m)
#  
n <- 5 # First n columns
A <- U[,1:n]
}

