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
library(penaltyLearning)
data(neuroblastomaProcessed, package="penaltyLearning")
library(doParallel)
registerDoParallel()

errors.per.model <- data.table(neuroblastomaProcessed$errors)
errors.per.model[, pid.chr := paste0(profile.id, ".", chromosome)]
setkey(errors.per.model, pid.chr)
set.seed(1)
n.train <- 500
i.train <- 1:n.train
fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
  feature.mat[i.train,], target.mat[i.train,],
  incorrect.labels.db=errors.per.model))
fit$plot

if(require(iregnet)){
  data("penalty.learning", package="iregnet")
  set.seed(1)
  is.test <- grepl("chr1:", rownames(penalty.learning$X.mat))
  pfit <- with(penalty.learning, IntervalRegressionCV(
    X.mat[!is.test,], y.mat[!is.test,]))
  print(pfit$plot)
  pred.log.lambda <- pfit$predict(penalty.learning$X.mat)
  residual <- targetIntervalResidual(penalty.learning$y.mat, pred.log.lambda)
  residual.tall <- data.table(is.test, residual)[, list(
    mean.residual=mean(residual),
    intervals=.N
    ), by=.(
         set=ifelse(is.test, "test", "train"),
         sign.residual=sign(residual))]
  residual.tall[, set.intervals := ifelse(
    set=="train", sum(!is.test), sum(is.test))]
  residual.tall[, percent.intervals := 100 * intervals / set.intervals]
  dcast(residual.tall, set ~ sign.residual, value.var="mean.residual")
  dcast(residual.tall, set ~ sign.residual, value.var="percent.intervals")
}

}
