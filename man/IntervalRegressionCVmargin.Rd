\name{IntervalRegressionCVmargin}
\alias{IntervalRegressionCVmargin}
\title{IntervalRegressionCVmargin}
\description{Use cross-validation to fit an L1-regularized linear interval
regression model by optimizing both margin and regularization
parameters. This function just calls IntervalRegressionCV with a
margin.vec parameter that is computed based on the finite target
interval limits.}
\usage{IntervalRegressionCVmargin(feature.mat, target.mat, 
    ...)}
\arguments{
  \item{feature.mat}{Numeric feature matrix, n observations x p features.}
  \item{target.mat}{Numeric target matrix, n observations x 2 limits.}
  \item{\dots}{Passed to IntervalRegressionCV.}
}

\value{Model fit list from IntervalRegressionCV.}

\author{Toby Dylan Hocking}




\examples{
if(interactive()){
  library(penaltyLearning)
  data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
  if(require(future)){
    plan(multiprocess)
  }
  set.seed(1)
  i.train <- 1:nrow(neuroblastomaProcessed$feature.mat)
  fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
    feature.mat[i.train,], target.mat[i.train,],
    verbose=1))
}
}
