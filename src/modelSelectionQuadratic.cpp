/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionQuadratic.h"
#include <math.h>
#include <stdio.h>

int modelSelectionQuadratic
(const double *L, const double *complexity_vec, int *n_models,
 int *K, double *b){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(L[i-1] <= L[i]){
      return ERROR_QUAD_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_QUAD_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0; // index of largest selected model.
  b[0] = INFINITY;
  K[0] = 0;
  double lambda, min_lambda, best_k;
  for(int t=1; t<N; t++){
    // In the pseudocode this is the start of Algorithm 2.
    min_lambda=INFINITY;
    for(int k=0; k<t; k++){
      lambda = (L[k]-L[t])/(complexity_vec[t]-complexity_vec[k]);
      if(lambda < min_lambda){
	min_lambda = lambda;
	M = k+1;
      }
    }
    // In the pseudocode this is the end of Algorithm 2.
    b[M] = lambda;
    K[M] = t;
  }
  *n_models = M; 
  return 0;
}
