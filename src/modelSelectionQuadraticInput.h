/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_QUAD_IN_LOSS_NOT_DECREASING 1
#define ERROR_QUAD_IN_COMPLEXITY_NOT_INCREASING 2

int modelSelectionQuadraticInput(const double*, const double*, int*,
		      int*, double*);
