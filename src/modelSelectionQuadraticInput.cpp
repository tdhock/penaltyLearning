/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionQuadraticInput.h"
#include <math.h>
#include <stdio.h>

int modelSelectionQuadraticInput
(const double *L, const double *complexity_vec, int *n_models,
 int *K, double *b){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(L[i-1] <= L[i]){
      return ERROR_QUAD_IN_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_QUAD_IN_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0, best_k; // index of largest selected model.
  b[0] = INFINITY;
  K[0] = 0;
  double lambda, min_lambda;
  for(int t=1; t<N; t++){
    // In the pseudocode this is the start of Algorithm 2.
    min_lambda=INFINITY;
    for(int k=0; k<t; k++){
      lambda = (L[k]-L[t])/(complexity_vec[t]-complexity_vec[k]);
      if(lambda < min_lambda){
	min_lambda = lambda;
	best_k = k;
      }
    }
    while(best_k < K[M])M--;
    M++;
    // In the pseudocode this is the end of Algorithm 2.
    b[M] = min_lambda;
    K[M] = t;
  }
  *n_models = M; 
  return 0;
}
