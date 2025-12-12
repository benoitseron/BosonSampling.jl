// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]

arma::cx_vec cxPermMinors(arma::cx_mat C) {
// Declarations
	int n = C.n_cols, m = C.n_rows; // m = n + 1
	if(m != n + 1) stop("Input matrix has incorrect dimensions");
    int i, j = 0, k;  
    bool s = true; 
    arma::cx_vec p(m), q(m), v(m);
    arma::cx_double t;
//
    arma::uvec d(n); d.fill(true);  // almost a boolean vector
//  
    if(n == 1) return arma::flipud(C); 
//
// g: initialize Gray code iterator
//
	arma::ivec g = arma::regspace< arma::ivec>(0,(n-1));	 
// Initialise delta sums v and partial products p 
//
    v = arma::sum(C,1)/2;
    q = arma::cumprod(v); t = v[m-1];
    p[m-1] = q[m-2];
    for(i = m-2; i > 0; i--){
        p[i] = t*q[i-1];
        t *= v[i];
    }
    p[0] = t;
	while(j < n-1){	
        if(d[j]) v -= C.col(j); else v += C.col(j);
        q = arma::cumprod(v); t = v[m-1];
        if(s){
            p[m-1] -= q[m-2];
            for(i = m-2; i > 0; i--){
                p[i] -= t*q[i-1];
                t *= v[i];
            }
            p[0] -= t;
        } else {
           p[m-1] += q[m-2];
           for(i = m-2; i > 0; i--){
                p[i] += t*q[i-1];
                t *= v[i];		
            }
            p[0] += t;
        }
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
