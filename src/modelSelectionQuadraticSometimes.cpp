/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionQuadraticSometimes.h"
#include <math.h>
#include <stdio.h>

int modelSelectionQuadraticSometimes
(const double *L, const double *complexity_vec, int *n_models,
 int *K, double *b, int *iterations){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(L[i-1] <= L[i]){
      return ERROR_QUAD_SOME_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_QUAD_SOME_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0, best_i; // index of largest selected model.
  b[0] = INFINITY;
  K[0] = 0;
  iterations[0]=0;
  double lambda, best_lambda;
  for(int t=1; t<N; t++){
    // In the pseudocode this is the start of Algorithm 2.
    best_lambda=INFINITY;
    best_i=0;
    iterations[t]=M;
    for(int i=0; i <= M; i++){
      int k=K[i];
      lambda = (L[k]-L[t])/(complexity_vec[t]-complexity_vec[k]);
      if(lambda < b[i] && i > best_i){
	best_lambda = lambda;
	best_i = i;
      }
    }
    M = best_i + 1;
    // In the pseudocode this is the end of Algorithm 2.
    b[M] = best_lambda;
    K[M] = t;
  }
  *n_models = M; 
  return 0;
}
