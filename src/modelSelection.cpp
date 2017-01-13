/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include <list>

int modelSelection
(const double *loss_vec, const double *complexity_vec, const int n_input,
 int *optimal_before_vec, double *lambda_vec){
  std::list<breakpoint> break_list;
  return 0;
}
