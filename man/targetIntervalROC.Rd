\name{targetIntervalROC}
\alias{targetIntervalROC}
\title{targetIntervalROC}
\description{Compute a ROC curve using a target interval matrix. WARNING: this
ROC curve is less detailed than the one you get from ROChange! Use
ROChange if possible.}
\usage{targetIntervalROC(target.mat, pred)}
\arguments{
  \item{target.mat}{n x 2 numeric matrix: target intervals of log(penalty) values that
yield minimal incorrect labels.}
  \item{pred}{numeric vector: predicted log(penalty) values.}
}

\value{list describing ROC curves, same as ROChange.}

\author{Toby Dylan Hocking}





