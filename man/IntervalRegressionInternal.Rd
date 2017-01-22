\name{IntervalRegressionInternal}
\alias{IntervalRegressionInternal}
\title{IntervalRegressionInternal}
\description{Solve the squared hinge loss interval regression problem for one
regularization parameter: w* = argmin_w L(w) + regularization *
||w||_1 where L(w) is the average squared hinge loss with respect
to the targets, and ||w||_1 is the L1-norm of the weight vector
(excluding the first element, which is the un-regularized
intercept or bias term). This function performs no scaling of
input features, and is meant for internal use only! To learn a
regression model, try IntervalRegressionCV or
IntervalRegressionUnregularized.}
\usage{IntervalRegressionInternal(features, targets, initial.param.vec, 
    regularization, threshold = 0.001, max.iterations = 1000, 
    weight.vec = NULL, Lipschitz = NULL, verbose = 2)}
\arguments{
  \item{features}{Scaled numeric feature matrix (problems x features). The first
column/feature should be all ones and will not be regularized.}
  \item{targets}{Numeric target matrix (problems x 2).}
  \item{initial.param.vec}{initial guess for weight vector (features).}
  \item{regularization}{Degree of L1-regularization.}
  \item{threshold}{When the stopping criterion gets below this threshold, the
algorithm stops and declares the solution as optimal.}
  \item{max.iterations}{Error if the algorithm has not found an optimal solution after
this many iterations.}
  \item{weight.vec}{A numeric vector of weights for each training example.}
  \item{Lipschitz}{A numeric scalar or NULL, which means to compute Lipschitz as the
mean of the squared L2-norms of the rows of the feature matrix.}
  \item{verbose}{Cat messages: for restarts and at the end if >= 1, and for every
iteration if >= 2.}
}

\value{Numeric vector of scaled weights w of the affine function f_w(X) =
X \%*\% w for a scaled feature matrix X with the first row entirely
ones.}

\author{Toby Dylan Hocking}





