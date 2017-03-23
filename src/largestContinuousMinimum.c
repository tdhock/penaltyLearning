/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "largestContinuousMinimum.h"
#include <math.h>

int largestContinuousMinimum
(int n_data, double *cost_vec, double *size_vec, int *index_vec){
  double largest_size = 0.0;
  double min_cost = -INFINITY;
  double cur_size = 0.0;
  double cur_cost = -INFINITY;
  int start;
  for(int i=0; i<n_data; i++){
    if(size_vec[i] <= 0.0){
      return ERROR_SIZES_MUST_BE_POSITIVE;
    }
    if(
       (i==0)//first data point.
       ||
       (cost_vec[i] != cur_cost)//different from prev cost.
       ){
      cur_size = size_vec[i];
      cur_cost = cost_vec[i];
      start = i;
    }else{
      cur_size += size_vec[i];
    }
    if(cur_cost <= min_cost && largest_size < cur_size){
      //store new min cost.
      min_cost = cur_cost;
      largest_size = cur_size;
      index_vec[0] = start;
      index_vec[1] = i;
    }
  }
  return 0;
}

