/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_RIGAILL_LOSS_NOT_DECREASING 1
#define ERROR_RIGAILL_COMPLEXITY_NOT_INCREASING 2

int modelSelectionRigaill(const double*, const double*, int*,
		      int*, double*);
