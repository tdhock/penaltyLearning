\name{largestContinuousMinimumR}
\alias{largestContinuousMinimumR}
\title{largestContinuousMinimumR}
\description{Find the run of minimum \code{cost} with the largest \code{size}.
This function uses a two pass R implementation,
and is meant for internal use.
Use \code{\link{targetIntervals}} for real data.}
\usage{largestContinuousMinimumR(cost, size)}
\arguments{
  \item{cost}{numeric vector of \code{cost} values.}
  \item{size}{numeric vector of interval \code{size} values.}
}

\value{Integer vector length 2 (start and end of target interval relative
to \code{cost} and \code{size}).}

\author{Toby Dylan Hocking}




\examples{

library(penaltyLearning)
data(neuroblastomaProcessed, envir=environment())
one.problem.error <-
  neuroblastomaProcessed$errors[profile.id=="4" & chromosome=="1"]
indices <- one.problem.error[, largestContinuousMinimumR(
  errors, max.log.lambda-min.log.lambda)]
one.problem.error[indices[["start"]]:indices[["end"]],]

}
