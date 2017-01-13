/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include <stdio.h>
#include <list>

int modelSelection
(const double *loss_vec, const double *complexity_vec, const int n_input,
 int *optimal_after_vec, double *lambda_vec){
  std::list<breakpoint> break_list;
  int lo=n_input-2, hi=n_input-1;
  double last_lambda = 0;
  while(0 <= lo){
    double lambda =
      (loss_vec[lo]-loss_vec[hi])/(complexity_vec[hi]-complexity_vec[lo]);
    //printf("lo,hi=%d,%d last=%f lambda=%f\n", lo+1, hi+1, last_lambda, lambda);
    if(last_lambda < lambda){
      //printf("push(%f,%d)\n", lambda, hi+1);
      break_list.emplace_back(lambda, lo);
      hi = lo;
      lo = hi - 1;
      last_lambda = lambda;
    }else{
      //printf("pop()\n");
      break_list.pop_back();
      if(break_list.size()){
	hi = break_list.back().optimal_after;
      }else{
	hi = n_input-1;
      }
      last_lambda = break_list.back().penalty;
    }
  }
  std::list<breakpoint>::iterator it;
  //printf("%d elements in list\n", break_list.size());
  for(it=break_list.begin(), lo=0; it != break_list.end(); it++, lo++){
    optimal_after_vec[lo] = it->optimal_after;
    lambda_vec[lo] = it->penalty;
  }
  return 0;
}
