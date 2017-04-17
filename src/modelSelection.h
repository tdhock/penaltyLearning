/* -*- compile-command: "R CMD INSTALL .." -*- */

#define ERROR_LOSS_NOT_DECREASING 1
#define ERROR_COMPLEXITY_NOT_INCREASING 2

int modelSelection(const double*, const double*, int,
		   int*, double*);

class breakpoint {
public:
  double penalty;
  int optimal_after;
  breakpoint(double pen, int after){
    penalty = pen;
    optimal_after = after;
  }
};

