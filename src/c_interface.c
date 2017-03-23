/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "largestContinuousMinimum.h"
#include <R.h>

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

