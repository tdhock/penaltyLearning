/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include <list>

int modelSelection
(double *loss_vec, double *complexity_vec, int n_input,
 int *optimal_before_vec, double *lambda_vec, int *n_output){
  std::list<breakpoint> break_list;
  return 0;
}
