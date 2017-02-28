\name{IntervalRegressionCV}
\alias{IntervalRegressionCV}
\title{IntervalRegressionCV}
\description{Use cross-validation to estimate the optimal regularization, by
picking the value that minimizes the number of incorrectly
predicted target intervals. K-fold cross-validation is
parallelized using the foreach package.}
\usage{IntervalRegressionCV(feature.mat, target.mat, n.folds = ifelse(nrow(feature.mat) < 
    10, 3L, 5L), fold.vec = sample(rep(1:n.folds, l = nrow(feature.mat))), 
    verbose = 0, min.observations = 10, reg.type = "1sd", 
    incorrect.labels.db = NULL)}
\arguments{
  \item{feature.mat}{Numeric feature matrix.}
  \item{target.mat}{Numeric target matrix.}
  \item{n.folds}{Number of cross-validation folds.}
  \item{fold.vec}{Integer vector of fold id numbers.}
  \item{verbose}{numeric: bigger numbers for more output.}
  \item{min.observations}{stop with an error if there are fewer than this many observations.}
  \item{reg.type}{Either "1sd", "min(mean)" or "mean(min)" which specifies how the
regularization parameter is chosen during the interval
cross-validation loop.}
  \item{incorrect.labels.db}{either NULL or a data.table, which specifies what to compute to
select the regularization parameter on the validation set. NULL
means to minimize the squared hinge loss, which measures how far
the predicted log(penalty) values are from the target
intervals. If a data.table is specified, it should have one key
which corresponds to the rownames of feature.mat, and columns
min.log.lambda, max.log.lambda, fp, fn, possible.fp, possible.fn;
these will be used with ROChange to compute the AUC for each
regularization parameter, and the maximimum will be selected. This
data.table can be computed via
labelError(modelSelection(...),...)$model.errors -- see
example(ROChange).}
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
fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
  feature.mat, target.mat,
  incorrect.labels.db=errors.per.model))
fit$plot

if(require(iregnet)){
  data("penalty.learning", package="iregnet")
  set.seed(1)
  is.test <- grepl("chr1:", rownames(penalty.learning$X.mat))
  pfit <- with(penalty.learning, IntervalRegressionCV(X.mat[!is.test,], y.mat[!is.test,]))
  print(pfit$plot)
  pred.log.lambda <- pfit$predict(penalty.learning$X.mat)
  residual <- targetIntervalResidual(penalty.learning$y.mat, pred.log.lambda)
  residual.tall <- data.table(is.test, residual)[, list(
    mean.residual=mean(residual),
    intervals=.N
    ), by=.(set=ifelse(is.test, "test", "train"), sign.residual=sign(residual))]
  residual.tall[, set.intervals := ifelse(set=="train", sum(!is.test), sum(is.test))]
  residual.tall[, percent.intervals := 100 * intervals / set.intervals]
  dcast(residual.tall, set ~ sign.residual, value.var="mean.residual")
  dcast(residual.tall, set ~ sign.residual, value.var="percent.intervals")
}

}
