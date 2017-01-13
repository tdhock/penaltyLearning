/* -*- compile-command: "R CMD INSTALL .." -*- */

int modelSelection(const double*, const double*, int,
		   int*, double*);

class breakpoint {
public:
  double penalty;
  int optimal_before;
  breakpoint(double pen, int i){
    penalty = pen;
    optimal_before = i;
  }
};

