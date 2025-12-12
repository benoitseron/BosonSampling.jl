// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

arma::cx_double cxPerm(arma::cx_mat A) {
// Declarations
	int n = A.n_cols,  m = A.n_rows; // needs to be square matrix 	
	if(m != n) stop("Input matrix has incorrect dimensions");
    int j = 0, k; 
	arma::cx_vec v = arma::sum(A,1)/2.;
	arma::cx_double p = arma::prod(v);
	bool s = true;
//
	arma::uvec d(n); d.fill(true); // almost a boolean vector
// 
  if(n == 1) return A(0,0); 
//
// g: initialize Gray code iterator
	arma::ivec g = arma::regspace< arma::ivec >(0,(n-1));
//
	while(j < n-1){
		if(d[j]) v -= A.col(j); else v += A.col(j);
		s ? p -= arma::prod(v) : p += arma::prod(v);
		d[j] = !d[j]; s = !s;   
// iterate Gray code: j is active index  		
		if( j > 0){
            k = j + 1; g[j] = g[k]; g[k] = k; j = 0;
		} else {
			j = g[1]; g[1] = 1;
		}  	
	}
    return 2.*p;
}

//END
