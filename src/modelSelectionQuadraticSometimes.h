/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_QUAD_SOME_LOSS_NOT_DECREASING 1
#define ERROR_QUAD_SOME_COMPLEXITY_NOT_INCREASING 2

int modelSelectionQuadraticSometimes
(const double*, const double*, int*,
 int*, double*, int*);
