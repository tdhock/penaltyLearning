/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelectionLinear.h"
#include <math.h>
#include <stdio.h>

// This is the c(t,i) function from the paper, which returns the
// penalty/lambda value for which f_t(penalty) =
// f_{K_{t-1,i}}(penalty), where f_k(penalty)=L_k+penalty*C_k (in the
// paper C[k]=k for all k from 1 to N).
double c(int t, int i, int *K, const double *L, const double *C){
  int k = K[i];
  return (L[k]-L[t])/(C[t]-C[k]);
}

// Return value is an error status code defined in
// modelSelectionLinear.h.
int modelSelectionLinear
// This is Algorithm 1 from "Linear time dynamic programming for the
// exact path of optimal models selected from a finite set." Argument
// names in this C code are consistent with the variable names /
// pseudocode in the paper. This file uses the .cpp extension so that
// it interfaces easily with R and the other compiled code in this
// package, but since it is just standard C you should be able to copy
// it and use in other contexts with a standard C compiler. All
// arguments below are pointers to memory which should be
// pre-allocated before calling this function.
(const double *L, //Input[N]: decreasing loss values.
 const double *complexity_vec,  //Input[N]:
 // increasing complexity values. In the paper this is fixed as
 // complexity_vec[t]=t+1 for all t from 0 to N-1 (actually t=1
 // to N in the 1-indexed paper/math but C code is 0-indexed), which
 // corresponds to using the L0 pseudo-norm as the model complexity
 // (number of non-zero elements). This is included as an option in
 // order to support other measures of model complexity.
 int *n_models, //Input/Output[1]:
 // number of loss/complexity values at input (N from math/paper),
 // number of selected models at output.
 int *K, //Output[N]: selected models,
 // i.e. at the end of the function, for all i from 0 to N-1, K[i] in
 // the code is the same as K_{N,i+1} in the math/paper.
 double *b, //Output[N]: breakpoints,
 // i.e. at the end of the function, for all i from 0 to N-1, b[i] in
 // the code is b_{N,i} in the math/paper. Note that the first
 // breakpoint b_{N,0}=INFINITY is stored in b[0], but the last
 // breakpoint b_{N,M_N}=0 is NOT.
 int *w //Output[N]: while loop iterations.
 // not required for computing the optimal models/breakpoints, but
 // stored in order to analyze the empirical time complexity of the
 // algorithm.
 ){
  int N = *n_models;
  for(int i=1; i < N; i++){
    if(L[i-1] <= L[i]){
      return ERROR_LINEAR_LOSS_NOT_DECREASING;
    }
    if(complexity_vec[i] <= complexity_vec[i-1]){ 
      return ERROR_LINEAR_COMPLEXITY_NOT_INCREASING;
    }
  }
  int M = 0; // index of largest selected model.
  b[0] = INFINITY;
  K[0] = 0;
  double lambda;
  for(int t=1; t<N; t++){
    // In the pseudocode this is the start of Algorithm 2.
    w[t] = 1;
    while
      ((lambda = c(t, M, K, L, complexity_vec)) >=
       b[M]){
      M--;
      w[t]++;
    }
    M++;
    // In the pseudocode this is the end of Algorithm 2.
    b[M] = lambda;
    K[M] = t;
  }
  *n_models = M; 
  return 0;
}
