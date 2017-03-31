\name{targetIntervalROC}
\alias{targetIntervalROC}
\title{targetIntervalROC}
\description{Compute a ROC curve using a target interval matrix. A prediction
less than the lower limit is considered a false positive (penalty
too small, too many changes), and a prediction greater than the
upper limit is a false negative (penalty too large, too few
changes). WARNING: this ROC curve is less detailed than the one
you get from ROChange! Use ROChange if possible.}
\usage{targetIntervalROC(target.mat, pred)}
\arguments{
  \item{target.mat}{n x 2 numeric matrix: target intervals of log(penalty) values that
yield minimal incorrect labels.}
  \item{pred}{numeric vector: predicted log(penalty) values.}
}

\value{list describing ROC curves, same as ROChange.}

\author{Toby Dylan Hocking}




\examples{

## TODO plot that compares ROChange and targetIntervalROC error
## curves.

}
