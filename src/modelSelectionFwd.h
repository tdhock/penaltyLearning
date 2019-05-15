/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_FWD_LOSS_NOT_DECREASING 1
#define ERROR_FWD_COMPLEXITY_NOT_INCREASING 2

int modelSelectionFwd(const double*, const double*, int*,
		      int*, double*, int*);
