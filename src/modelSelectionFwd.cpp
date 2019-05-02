/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionFwd.h"
#include <stdio.h>

int modelSelectionFwd
(const double *loss_vec, const double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec, int *loop_eval_vec){
  for(int i=1; i < *n_models; i++){
    if(loss_vec[i-1] <= loss_vec[i]){
      return ERROR_FWD_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_FWD_COMPLEXITY_NOT_INCREASING;
    }
  }
  return 0;
}
