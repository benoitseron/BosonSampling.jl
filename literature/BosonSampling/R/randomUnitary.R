#' Random unitary
#'
#' Generates a random unitary matrix (square)
#' @param size Dimension of matrix 
#' @keywords complex
#' @export
#' @examples
#' m = 25
#' set.seed(7)
#' U = randomUnitary(m)
#' #  
#' n = 5 # First n columns
#' A = U[,1:n]

randomUnitary = function(size){
	B = qr(matrix(complex(real = stats::rnorm(size^2), imaginary = stats::rnorm(size^2)),size))
	R_Diag = sign(Re(diag(qr.R(B))))
  qr.Q(B)%*%diag(R_Diag)
}

# ## 
# k = 25
# A = matrix(complex(real = rnorm(k^2), imaginary = rnorm(k^2)),nrow = k,byrow = F)
# B = as.matrix(A[,-k])
# c(cxPerm(A),sum(cxPermMinorsV2(B)*A[,k])) 