/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionQuadraticOutput.h"
#include <math.h>
#include <stdio.h>

int modelSelectionQuadraticOutput
(const double *L, const double *complexity_vec, int *n_models,
 int *K, double *b){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(L[i-1] <= L[i]){
      return ERROR_QUAD_OUT_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_QUAD_OUT_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0, best_k=0; // index of largest selected model.
  b[0] = INFINITY;
  K[0] = 0;
  double lambda, best_lambda;
  for(int t=1; t<N; t++){
    // In the pseudocode this is the start of Algorithm 2.
    best_lambda=INFINITY;
    for(int i=0; i <= M; i++){
      int k=K[i];
      lambda = (L[k]-L[t])/(complexity_vec[t]-complexity_vec[k]);
      if(lambda < best_lambda && k > best_k){
	best_lambda = lambda;
	best_k = k;
      }
    }
    // In the pseudocode this is the end of Algorithm 2.
    b[M] = best_lambda;
    K[M] = t;
  }
  *n_models = M; 
  return 0;
}
