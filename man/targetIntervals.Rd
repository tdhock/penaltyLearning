\name{targetIntervals}
\alias{targetIntervals}
\title{Compute target intervals}
\description{Compute target intervals of log(penalty) values that result in
predicted changepoint \code{models} with minimum incorrect labels.
Use this function after \code{\link{labelError}}, and before IntervalRegression*.}
\usage{targetIntervals(models, 
    problem.vars)}
\arguments{
  \item{models}{data.table with columns errors, min.log.lambda, max.log.lambda,
typically labelError()$model.errors.}
  \item{problem.vars}{character: column names used to identify data set / segmentation
problem.}
}

\value{data.table with columns \code{problem.vars}, one row for each
segmentation problem. The "min.log.lambda", and "max.log.lambda"
columns give the largest interval of log(penalty) values which
results in the minimum incorrect labels for that problem. This can
be used to create the target.mat parameter of the
IntervalRegression* functions.}

\author{Toby Dylan Hocking}




\examples{

data.table::setDTthreads(1)

library(penaltyLearning)
data(neuroblastomaProcessed, envir=environment())
targets.dt <- targetIntervals(
  neuroblastomaProcessed$errors,
  problem.vars=c("profile.id", "chromosome"))

}
