/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include <R.h>

extern "C" {

  void modelSelection_interface
  (double *loss, double *complexity, int *n_models,
   int *before, double *lambda
   ){
    int status = modelSelection(loss, complexity, *n_models, before, lambda);
    if(status == ERROR_LOSS_NOT_DECREASING){
      error("loss not decreasing");
    }
    if(status == ERROR_COMPLEXITY_NOT_INCREASING){
      error("complexity not increasing");
    }
    if(status != 0){
      error("error code %d", status);
    }
  }
  
}
    
