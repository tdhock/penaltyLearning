/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_QUAD_OUT_LOSS_NOT_DECREASING 1
#define ERROR_QUAD_OUT_COMPLEXITY_NOT_INCREASING 2

int modelSelectionQuadraticOutput
(const double*, const double*, int*,
 int*, double*);
