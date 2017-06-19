\name{largestContinuousMinimumC}
\alias{largestContinuousMinimumC}
\title{largestContinuousMinimumC}
\description{Find the run of minimum cost with the largest size.
This function use a linear time C implementation,
and is meant for internal use.
Use targetIntervals for real data.}
\usage{largestContinuousMinimumC(cost, size)}
\arguments{
  \item{cost}{numeric vector of cost values.}
  \item{size}{numeric vector of interval size values.}
}

\value{Integer vector length 2 (start and end of target interval relative
to cost and size).}

\author{Toby Dylan Hocking}




\examples{

library(penaltyLearning)
data(neuroblastomaProcessed, envir=environment())
one.problem.error <-
  neuroblastomaProcessed$errors[profile.id=="4" & chromosome=="1"]
indices <- one.problem.error[, largestContinuousMinimumC(
  errors, max.log.lambda-min.log.lambda)]
one.problem.error[indices[["start"]]:indices[["end"]],]

}
