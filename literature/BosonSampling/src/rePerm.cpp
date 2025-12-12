// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

double rePerm(arma::mat B) {
// Declarations
	int n = B.n_cols, m = B.n_rows; // needs to be square matrix 
	if(m != n) stop("Input matrix has incorrect dimensions");
	int j = 0, k;
	arma::vec v = arma::sum(B,1)/2.;
	double p = arma::prod(v);
    bool s = true;
//	
    arma::uvec d(n); d.fill(true);
//
	if(n == 1) return B(0,0);
//
// g: initialize Gray code iterator
	arma::ivec g = arma::regspace< arma::ivec >(0,(n-1));
//	
	while(j < n-1){	 	
		if(d[j]) v -= B.col(j); else v += B.col(j);
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

// END
