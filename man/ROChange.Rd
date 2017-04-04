\name{ROChange}
\alias{ROChange}
\title{ROC curve for changepoints}
\description{Compute a Receiver Operating Characteristic curve for a penalty
function.}
\usage{ROChange(models, predictions, problem.vars = character())}
\arguments{
  \item{models}{data.frame describing the number of incorrect labels as a function
of log(lambda), with columns min.log.lambda, max.log.lambda, fp,
fn, possible.fp, possible.fn, etc. This can be computed via
labelError(modelSelection(...), ...)$model.errors -- see examples.}
  \item{predictions}{data.frame with a column named pred.log.lambda, the predicted
log(lambda) value for each segmentation problem.}
  \item{problem.vars}{character: column names used to identify data set / segmentation
problem. }
}

\value{list of results describing ROC curve: roc is a data.table with one
row for each point on the ROC curve; thresholds is the two rows of
roc which correspond to the predicted and minimal error
thresholds; auc.polygon is a data.table with one row for each
vertex of the polygon used to compute AUC; auc is the numeric Area
Under the ROC curve, actually computed via geometry::polyarea as
the area inside the auc.polygon.}

\author{Toby Dylan Hocking}




\examples{

library(penaltyLearning)
data(neuroblastomaProcessed, envir=environment())
## Get incorrect labels data for one profile.
pid <- 11
pro.errors <- neuroblastomaProcessed$errors[profile.id==pid,]
## Get the feature that corresponds to the BIC penalty = log(n),
## meaning log(penalty) = log(log(n)).
chr.vec <- paste(c(1:4, 11, 17))
pid.names <- paste0(pid, ".", chr.vec)
BIC.feature <- neuroblastomaProcessed$feature.mat[pid.names, "log2.n"]
pred <- data.table(pred.log.lambda=BIC.feature, chromosome=chr.vec)
result <- ROChange(pro.errors, pred, "chromosome")
library(ggplot2)
## Plot the ROC curves.
ggplot()+
  geom_path(aes(FPR, TPR), data=result$roc)+
  geom_point(aes(FPR, TPR, color=threshold), data=result$thresholds, shape=1)

## Plot the number of incorrect labels as a function of threshold.
ggplot()+
  geom_segment(aes(
    min.thresh, errors,
    xend=max.thresh, yend=errors),
    data=result$roc)+
  geom_point(aes((min.thresh+max.thresh)/2, errors, color=threshold),
             data=result$thresholds,
             shape=1)+
  xlab("log(penalty) constant added to BIC penalty")

}
