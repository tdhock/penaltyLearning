/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionRigaill.h"
#include <math.h>
#include <stdio.h>

int modelSelectionRigaill
(const double *loss_vec, const double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(loss_vec[i-1] <= loss_vec[i]){
      return ERROR_RIGAILL_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_RIGAILL_COMPLEXITY_NOT_INCREASING;
    }
  }
  int kc = N-1, model_i=0;
  selected_model_vec[0]=N-1;
  selected_penalty_vec[0]=0;
  while(0 < kc){
    double next_lambda = INFINITY;
    int next_K = 0;
    for(int k=0; k<kc; k++){
      double hit_time =
	(loss_vec[kc] - loss_vec[k])/(complexity_vec[k]-complexity_vec[kc]);
      if(hit_time < next_lambda){
	next_lambda = hit_time;
	next_K = k;
      }
    }
    kc = next_K;
    model_i++;
    selected_model_vec[model_i]=kc;
    selected_penalty_vec[model_i]=next_lambda;
  }
  *n_models = model_i; 
  return 0;
}
