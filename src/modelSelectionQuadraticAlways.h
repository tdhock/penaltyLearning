/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_QUAD_ALWAYS_LOSS_NOT_DECREASING 1
#define ERROR_QUAD_ALWAYS_COMPLEXITY_NOT_INCREASING 2

int modelSelectionQuadraticAlways(const double*, const double*, int*,
				  int*, double*, int*);
