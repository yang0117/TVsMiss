#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

SEXP loglikelihood(SEXP beta_est_, SEXP data_matrix_){

    int n = nrows(data_matrix_);
    int p = ncols(data_matrix_);
    int lbeta = length(beta_est_);
    double mean_value=0.0;

    if(lbeta != p)
        error("x and beta_est does not match(loglikelihood)");
    beta_est_ = coerceVector(beta_est_, REALSXP);
    data_matrix_ = coerceVector(data_matrix_, REALSXP);

    double *beta_est = REAL(beta_est_);
    double *data_matrix = REAL(data_matrix_);


    int vec_length = n*(n-1)/2;
    double res[vec_length];
    for(int i=0; i < vec_length; i++){
      res[i] = 0;
    }


    int store_idx = 0;
    for(int i=0;i < n;i++ ){
      for(int j=i+1;j < n;j++){

        for(int k=0;k<p;k++){
          if (k == 0){
            res[store_idx] = beta_est[k]*0;
          }else{
            res[store_idx] += beta_est[k]*(data_matrix[i+k*n] - data_matrix[j+k*n]);
          }
        }
        res[store_idx] = -log(1+exp(-res[store_idx]*(data_matrix[i]-data_matrix[j])));
        store_idx++;
      }

    }


    //Rprintf("(%d)\n",vec_length);
    //Rprintf("(%d)\n",n);
    for(int i=0;i<vec_length;i++){
      mean_value += res[i];
    }
    mean_value = mean_value/vec_length;


    SEXP result;
    result = PROTECT(allocVector(REALSXP, 1));
    REAL(result)[0] = mean_value;
    UNPROTECT(1);
    return result;
}
