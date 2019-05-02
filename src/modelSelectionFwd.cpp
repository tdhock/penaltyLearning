/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionFwd.h"
#include <math.h>
#include <stdio.h>

double c(int t, int i, int *K, const double *L, const double *C){
  int k = K[i];
  return (L[k]-L[t])/(C[t]-C[k]);
}

int modelSelectionFwd
(const double *loss_vec, const double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec, int *loop_eval_vec){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(loss_vec[i-1] <= loss_vec[i]){
      return ERROR_FWD_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_FWD_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0; // index of largest selected model.
  selected_penalty_vec[0] = INFINITY;
  selected_model_vec[0] = 0;
  double lambda;
  for(int t=1; t<N; t++){
    loop_eval_vec[t] = 1;
    //printf("t=%d\n", t);
    while
      ((lambda = c(t, M, selected_model_vec, loss_vec, complexity_vec)) >
       selected_penalty_vec[M]){
      //printf("lambda=%f M=%d\n", lambda, M);
      M--;
      loop_eval_vec[t]++;
    }
    M++;
    selected_penalty_vec[M] = lambda;
    selected_model_vec[M] = t;
  }
  *n_models = M; 
  return 0;
}
