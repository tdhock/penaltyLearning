\name{targetIntervalResidual}
\alias{targetIntervalResidual}
\title{targetIntervalResidual}
\description{Compute residual of predicted penalties with respect to target
intervals.}
\usage{targetIntervalResidual(target.mat, pred)}
\arguments{
  \item{target.mat}{n x 2 numeric matrix: target intervals of log(penalty) values that
yield minimal incorrect labels.}
  \item{pred}{numeric vector: predicted log(penalty) values.}
}

\value{numeric vector of n residuals. Predictions that are too high
(above target.mat[,2]) get positive residuals (too few
changepoints), and predictions that are too low (below
target.mat[,1]) get negative residuals.}

\author{Toby Dylan Hocking}





