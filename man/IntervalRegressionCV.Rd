\name{IntervalRegressionCV}
\alias{IntervalRegressionCV}
\title{IntervalRegressionCV}
\description{Use cross-validation to estimate the optimal regularization, by
picking the value that minimizes the number of incorrectly
predicted target intervals. K-fold cross-validation is
parallelized using the foreach package.}
\usage{IntervalRegressionCV(feature.mat, target.mat, n.folds = ifelse(nrow(feature.mat) < 
    10, 3L, 5L), fold.vec = sample(rep(1:n.folds, l = nrow(feature.mat))), 
    verbose = 0, min.observations = 10, reg.type = "min(mean)", 
    incorrect.labels.db = NULL, initial.regularization = 0.001)}
\arguments{
  \item{feature.mat}{Numeric feature matrix.}
  \item{target.mat}{Numeric target matrix.}
  \item{n.folds}{Number of cross-validation folds.}
  \item{fold.vec}{Integer vector of fold id numbers.}
  \item{verbose}{numeric: bigger numbers for more output.}
  \item{min.observations}{stop with an error if there are fewer than this many observations.}
  \item{reg.type}{Either "1sd", "min(mean)" or "mean(min)" which specifies how the
regularization parameter is chosen during the internal
cross-validation loop. min(mean): first take the mean of the K-CV
error functions, then minimize it (this is the default since it
tends to yield the least test error). 1sd: take the least complex
model which is within one standard deviation of that minimum (this
model is typical a bit less accurate, but much less complex, so
better if you want to interpret the coefficients). mean(min): take
the min of each K-CV error function, and then take their mean.}
  \item{incorrect.labels.db}{either NULL or a data.table, which specifies the error function to
compute for selecting the regularization parameter on the
validation set. NULL means to minimize the squared hinge loss,
which measures how far the predicted log(penalty) values are from
the target intervals. If a data.table is specified, its first key
should correspond to the rownames of feature.mat, and columns
min.log.lambda, max.log.lambda, fp, fn, possible.fp, possible.fn;
these will be used with ROChange to compute the AUC for each
regularization parameter, and the maximimum will be selected (in
the plot this is negative.auc, which is minimized). This
data.table can be computed via
labelError(modelSelection(...),...)$model.errors -- see
example(ROChange).}
  \item{initial.regularization}{Passed to IntervalRegressionRegularized.}
}



\author{Toby Dylan Hocking}




\examples{

if(interactive()){
  library(penaltyLearning)
  data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
  if(require(doParallel)){
    registerDoParallel()
  }
  set.seed(1)
  i.train <- 1:200
  fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
    feature.mat[i.train,], target.mat[i.train,]))
  ## When only features and target matrices are specified for
  ## training, the squared hinge loss is used as the metric to
  ## minimize on the validation set.
  plot(fit)
  ## Create an incorrect labels data.table (first key is same as
  ## rownames of feature.mat and target.mat).
  errors.per.model <- data.table(neuroblastomaProcessed$errors)
  errors.per.model[, pid.chr := paste0(profile.id, ".", chromosome)]
  setkey(errors.per.model, pid.chr)
  set.seed(1)
  fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
    feature.mat[i.train,], target.mat[i.train,],
    ## The incorrect.labels.db argument is optional, but can be used if
    ## you want to use AUC as the CV model selection criterion.
    incorrect.labels.db=errors.per.model))
  plot(fit)
}

}
