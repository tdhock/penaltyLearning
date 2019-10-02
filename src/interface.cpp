/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include "modelSelectionLinear.h"
#include "modelSelectionQuadraticAlways.h"
#include "modelSelectionQuadraticSometimes.h"
#include "modelSelectionRigaill.h"
#include "largestContinuousMinimum.h"
#include <R.h>
#include <R_ext/Rdynload.h>

void largestContinuousMinimum_interface
(int *n_data, double *cost_vec, double *size_vec, int *index_vec){
  int status = largestContinuousMinimum
    (*n_data, cost_vec, size_vec, index_vec);
  if(status==ERROR_SIZES_MUST_BE_POSITIVE){
    error("sizes must be positive");
  }
  if(status != 0){
    error("error code %d", status);
  }
}

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
   
void modelSelectionLinear_interface
(double *loss_vec, double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec,
 int *loop_eval_vec
 ){
  int status = modelSelectionLinear
    (loss_vec, complexity_vec, n_models,
     selected_model_vec, selected_penalty_vec, loop_eval_vec);
  if(status == ERROR_LINEAR_LOSS_NOT_DECREASING){
    error("loss not decreasing");
  }
  if(status == ERROR_LINEAR_COMPLEXITY_NOT_INCREASING){
    error("complexity not increasing");
  }
  if(status != 0){
    error("error code %d", status);
  }
}
  
void modelSelectionQuadraticAlways_interface
(double *loss_vec, double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec, int *iterations
 ){
  int status = modelSelectionQuadraticAlways
    (loss_vec, complexity_vec, n_models,
     selected_model_vec, selected_penalty_vec, iterations);
  if(status == ERROR_QUAD_ALWAYS_LOSS_NOT_DECREASING){
    error("loss not decreasing");
  }
  if(status == ERROR_QUAD_ALWAYS_COMPLEXITY_NOT_INCREASING){
    error("complexity not increasing");
  }
  if(status != 0){
    error("error code %d", status);
  }
}

void modelSelectionQuadraticSometimes_interface
(double *loss_vec, double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec, int *iterations
 ){
  int status = modelSelectionQuadraticSometimes
    (loss_vec, complexity_vec, n_models,
     selected_model_vec, selected_penalty_vec, iterations);
  if(status == ERROR_QUAD_SOME_LOSS_NOT_DECREASING){
    error("loss not decreasing");
  }
  if(status == ERROR_QUAD_SOME_COMPLEXITY_NOT_INCREASING){
    error("complexity not increasing");
  }
  if(status != 0){
    error("error code %d", status);
  }
}

void modelSelectionRigaill_interface
(double *loss_vec, double *complexity_vec, int *n_models,
 int *selected_model_vec, double *selected_penalty_vec, int *iterations
 ){
  int status = modelSelectionRigaill
    (loss_vec, complexity_vec, n_models,
     selected_model_vec, selected_penalty_vec, iterations);
  if(status == ERROR_RIGAILL_LOSS_NOT_DECREASING){
    error("loss not decreasing");
  }
  if(status == ERROR_RIGAILL_COMPLEXITY_NOT_INCREASING){
    error("complexity not increasing");
  }
  if(status != 0){
    error("error code %d", status);
  }
}
 
R_CMethodDef cMethods[] = {
  {"modelSelectionQuadraticAlways_interface",
   (DL_FUNC) &modelSelectionQuadraticAlways_interface, 6
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"modelSelectionQuadraticSometimes_interface",
   (DL_FUNC) &modelSelectionQuadraticSometimes_interface, 6
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"modelSelectionRigaill_interface",
   (DL_FUNC) &modelSelectionRigaill_interface, 6
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"modelSelectionLinear_interface",
   (DL_FUNC) &modelSelectionLinear_interface, 6
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"modelSelection_interface",
   (DL_FUNC) &modelSelection_interface, 5
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"largestContinuousMinimum_interface",
   (DL_FUNC) &largestContinuousMinimum_interface, 4
   //,{INTSXP, REALSXP, REALSXP, INTSXP}
  },
  {NULL, NULL, 0}
};

extern "C" {
  void R_init_penaltyLearning(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    //R_useDynamicSymbols call says the DLL is not to be searched for
    //entry points specified by character strings so .C etc calls will
    //only find registered symbols.
    R_useDynamicSymbols(info, FALSE);
  }
}
