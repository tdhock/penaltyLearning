/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionQuadraticAlways.h"
#include <math.h>
#include <stdio.h>

int modelSelectionQuadraticAlways
(const double *L, const double *complexity_vec, int *n_models,
 int *K, double *b, int *iterations){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(L[i-1] <= L[i]){
      return ERROR_QUAD_ALWAYS_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_QUAD_ALWAYS_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0, best_k; // index of largest selected model.
  b[0] = INFINITY;
  K[0] = 0;
  iterations[0]=0;
  double lambda, best_lambda;
  for(int t=1; t<N; t++){
    // In the pseudocode this is the start of Algorithm 2.
    best_lambda=INFINITY;
    best_k=0;
    iterations[t]=t;
    for(int k=0; k<t; k++){
      lambda = (L[k]-L[t])/(complexity_vec[t]-complexity_vec[k]);
      if(lambda < best_lambda){
	best_lambda = lambda;
	best_k = k;
      }
    }
    while(best_k < K[M])M--;
    M++;
    // In the pseudocode this is the end of Algorithm 2.
    b[M] = best_lambda;
    K[M] = t;
  }
  *n_models = M; 
  return 0;
}
