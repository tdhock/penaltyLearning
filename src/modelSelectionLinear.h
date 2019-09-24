/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_LINEAR_LOSS_NOT_DECREASING 1
#define ERROR_LINEAR_COMPLEXITY_NOT_INCREASING 2

int modelSelectionLinear
(const double*, const double*, int*,
 int*, double*, int*);
