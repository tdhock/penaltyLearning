\name{targetIntervalROC}
\alias{targetIntervalROC}
\title{targetIntervalROC}
\description{Compute a ROC curve using a target interval matrix. A prediction
less than the lower limit is considered a false positive (penalty
too small, too many changes), and a prediction greater than the
upper limit is a false negative (penalty too large, too few
changes). WARNING: this ROC curve is less detailed than the one
you get from \code{\link{ROChange}}! Use \code{\link{ROChange}} if possible.}
\usage{targetIntervalROC(target.mat, pred)}
\arguments{
  \item{target.mat}{n x 2 numeric matrix: target intervals of log(penalty) values that
yield minimal incorrect labels.}
  \item{pred}{numeric vector: predicted log(penalty) values.}
}

\value{list describing ROC curves, same as \code{\link{ROChange}}.}

\author{Toby Dylan Hocking}




\examples{

library(penaltyLearning)
library(data.table)
data(neuroblastomaProcessed, envir=environment())

pid.vec <- c("1", "4")
chr <- 2
incorrect.labels <-
  neuroblastomaProcessed$errors[profile.id\%in\%pid.vec & chromosome==chr]
pid.chr <- paste0(pid.vec, ".", chr)
target.mat <- neuroblastomaProcessed$target.mat[pid.chr, , drop=FALSE]
pred.dt <- data.table(profile.id=pid.vec, pred.log.lambda=1.5)
roc.list <- list(
  labels=ROChange(incorrect.labels, pred.dt, "profile.id"),
  targets=targetIntervalROC(target.mat, pred.dt$pred.log.lambda))

err <- data.table(incorrect=names(roc.list))[, {
  roc.list[[incorrect]]$roc
}, by=incorrect]
library(ggplot2)
ggplot()+
  ggtitle("incorrect targets is an approximation of incorrect labels")+
  scale_size_manual(values=c(labels=2, targets=1))+
  geom_segment(aes(
    min.thresh, errors,
    color=incorrect,
    size=incorrect,
    xend=max.thresh, yend=errors),
               data=err)

}
