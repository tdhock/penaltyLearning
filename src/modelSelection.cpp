/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include <stdio.h>
#include <list>

int modelSelection
(const double *loss_vec, const double *complexity_vec, const int n_input,
 int *optimal_before_vec, double *lambda_vec){
  std::list<breakpoint> break_list;
  int lo=n_input-2, hi=n_input-1;
  double last_lambda = 0;
  while(0 <= lo){
    double lambda =
      (loss_vec[lo]-loss_vec[hi])/(complexity_vec[hi]-complexity_vec[lo]);
    //printf("lo,hi=%d,%d\n", lo, hi);
    if(last_lambda < lambda){
      //printf("push(%f,%d)\n", lambda, hi);
      break_list.emplace_back(lambda, hi);
      hi = lo;
      lo = hi - 1;
      last_lambda = lambda;
    }else{
      //printf("pop()\n");
      break_list.pop_back();
      hi++;
      last_lambda = break_list.back().penalty;
    }
  }
  std::list<breakpoint>::iterator it;
  for(it=break_list.begin(), lo=0; it != break_list.end(); it++, lo++){
    optimal_before_vec[lo] = it->optimal_before;
    lambda_vec[lo] = it->penalty;
  }
  return 0;
}
