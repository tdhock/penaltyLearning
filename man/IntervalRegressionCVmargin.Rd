\name{IntervalRegressionCVmargin}
\alias{IntervalRegressionCVmargin}
\title{IntervalRegressionCVmargin}
\description{Use cross-validation to fit an L1-regularized linear interval
regression model by optimizing both margin and regularization
parameters. This function just calls \code{\link{IntervalRegressionCV}} with a
margin.vec parameter that is computed based on the finite target
interval limits. If default parameters are used, this function
should be about 10 times slower than \code{\link{IntervalRegressionCV}}
(since this function computes n.margin=10 models
per regularization parameter whereas \code{\link{IntervalRegressionCV}}
only computes one).
On large (N > 1000 rows) data sets,
this function should yield a model which is a little
more accurate than \code{\link{IntervalRegressionCV}}
(since the margin parameter is optimized).}
\usage{IntervalRegressionCVmargin(feature.mat, 
    target.mat, log10.diff = 2, 
    n.margin = 10L, ...)}
\arguments{
  \item{feature.mat}{Numeric feature matrix, n observations x p features.}
  \item{target.mat}{Numeric target matrix, n observations x 2 limits.}
  \item{log10.diff}{Numeric scalar: factors of 10 below the largest finite limit
difference to use as a minimum margin value (difference on the
log10 scale which is used to generate margin parameters). Bigger
values mean a grid of margin parameters with a larger range. For
example if the largest finite limit in \code{target.mat} is 26 and the
smallest finite limit is -4 then the largest limit difference is
30, which will be used as the maximum margin parameter. If
\code{log10.diff} is the default of 2 then that means the smallest margin
parameter will be 0.3 (two factors of 10 smaller than 30).}
  \item{n.margin}{Integer scalar: number of margin parameters, by default 10.}
  \item{\dots}{Passed to \code{\link{IntervalRegressionCV}}.}
}

\value{Model fit list from \code{\link{IntervalRegressionCV}}.}

\author{Toby Dylan Hocking}




\examples{
if(interactive()){
  library(penaltyLearning)
  data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
  if(require(future)){
    plan(multiprocess)
  }
  set.seed(1)
  fit <- with(neuroblastomaProcessed, IntervalRegressionCVmargin(
    feature.mat, target.mat, verbose=1))
  plot(fit)
  print(fit$plot.heatmap)
}
}
