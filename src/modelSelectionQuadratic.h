/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_QUAD_LOSS_NOT_DECREASING 1
#define ERROR_QUAD_COMPLEXITY_NOT_INCREASING 2

int modelSelectionQuadratic(const double*, const double*, int*,
		      int*, double*);
