\name{IntervalRegressionCV}
\alias{IntervalRegressionCV}
\title{IntervalRegressionCV}
\description{Use cross-validation to fit an L1-regularized linear interval
regression model by optimizing margin and regularization
parameters (max AUC or min squared hinge loss). K-fold
cross-validation is parallelized using the future package. This
function repeatedly calls IntervalRegressionRegularized, and by
default assumes that margin=1 (to optimize the margin, specify it
manually, or use IntervalRegressionCVmargin).}
\usage{IntervalRegressionCV(feature.mat, target.mat, n.folds = ifelse(nrow(feature.mat) < 
    10, 3L, 5L), fold.vec = sample(rep(1:n.folds, l = nrow(feature.mat))), 
    verbose = 0, min.observations = 10, reg.type = "min", 
    incorrect.labels.db = NULL, initial.regularization = 0.001, 
    margin.vec = 1)}
\arguments{
  \item{feature.mat}{Numeric feature matrix, n observations x p features.}
  \item{target.mat}{Numeric target matrix, n observations x 2 limits.}
  \item{n.folds}{Number of cross-validation folds.}
  \item{fold.vec}{Integer vector of fold id numbers.}
  \item{verbose}{numeric: 0 for silent, bigger numbers (1 or 2) for more output.}
  \item{min.observations}{stop with an error if there are fewer than this many observations.}
  \item{reg.type}{Either "1sd" or "min" which specifies how the regularization
parameter is chosen during the internal cross-validation
loop. min: first take the mean of the K-CV error functions, then
minimize it (this is the default since it tends to yield the least
test error). 1sd: take the most regularized model with the same
margin which is within one standard deviation of that minimum 
(this model is typically a bit less accurate, but much less
complex, so better if you want to interpret the coefficients).}
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
  \item{margin.vec}{numeric vector of margin size hyper-parameters.}
}

\value{List representing regularized linear model.}

\author{Toby Dylan Hocking}




\examples{

if(interactive()){
  library(penaltyLearning)
  data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
  if(require(future)){
    plan(multiprocess)
  }
  set.seed(1)
  i.train <- 1:100
  fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
    feature.mat[i.train,], target.mat[i.train,],
    verbose=1))
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
