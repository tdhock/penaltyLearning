/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "largestContinuousMinimum.h"
#include <math.h>
#include <stdio.h>

int largestContinuousMinimum
(int n_data, double *cost_vec, double *size_vec, int *index_vec){
  double largest_size = 0.0;
  double min_cost = INFINITY;
  double cur_size = 0.0;
  double cur_cost = INFINITY;
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
    if(
       (cur_cost == min_cost && largest_size <= cur_size)//same cost but bigger
       ||
       (cur_cost < min_cost)//lower cost
       ){
      //store new min cost.
      //printf("cur_cost=%f min_cost=%f largest_size=%f cur_size=%f i=%d start=%d\n", cur_cost, min_cost, largest_size, cur_size, i, start);
      min_cost = cur_cost;
      largest_size = cur_size;
      index_vec[0] = start;
      index_vec[1] = i;
    }
  }
  //special case for 0 cost on both ends.
  if(min_cost == cost_vec[0]){
    index_vec[0] = 0;
  }
  return 0;
}

