\name{IntervalRegressionRegularized}
\alias{IntervalRegressionRegularized}
\title{IntervalRegressionRegularized}
\description{Repeatedly use \code{\link{IntervalRegressionInternal}} to solve interval
regression problems for a path of regularization parameters. This
function does not perform automatic selection of the
regularization parameter; instead, it returns regression models
for a range of regularization parameters, and it is up to you to
select which one to use. For automatic regularization parameter
selection, use \code{\link{IntervalRegressionCV}}.}
\usage{IntervalRegressionRegularized(feature.mat, target.mat, 
    initial.regularization = 0.001, factor.regularization = 1.2, 
    verbose = 0, margin = 1, ...)}
\arguments{
  \item{feature.mat}{Numeric feature matrix.}
  \item{target.mat}{Numeric target matrix.}
  \item{initial.regularization}{Initial regularization parameter.}
  \item{factor.regularization}{Increase regularization by this factor after finding an optimal
solution. Or NULL to compute just one model
(\code{initial.regularization}).}
  \item{verbose}{Print messages if >= 1.}
  \item{margin}{Non-negative \code{margin} size parameter, default 1.}
  \item{\dots}{Other parameters to pass to \code{\link{IntervalRegressionInternal}}.}
}

\value{List representing fit model. You can do
fit$predict(feature.matrix) to get a matrix of predicted log
penalty values. The param.mat is the n.features * n.regularization
numeric matrix of optimal coefficients (on the original scale).}

\author{Toby Dylan Hocking}




\examples{

if(interactive()){
  library(penaltyLearning)
  data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
  i.train <- 1:500
  fit <- with(neuroblastomaProcessed, IntervalRegressionRegularized(
    feature.mat[i.train,], target.mat[i.train,]))
  plot(fit)
}

}
