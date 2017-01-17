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
data(neuroblastoma, package="neuroblastoma", envir=environment())
pid <- 81
pro <- subset(neuroblastoma$profiles, profile.id==pid)
ann <- subset(neuroblastoma$annotations, profile.id==pid)
max.segments <- 20
segs.list <- list()
selection.list <- list()
for(chr in unique(ann$chromosome)){
  pro.chr <- subset(pro, chromosome==chr)
  fit <- Segmentor3IsBack::Segmentor(
    pro.chr$logratio, model=2, Kmax=max.segments)
  model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
  selection.df <- modelSelection(model.df, complexity="n.segments")
  selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
    (pro.chr$position[data.before.change]+
     pro.chr$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro.chr$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro.chr$position))
    segs.list[[paste(chr, n.segments)]] <- data.table(
      chromosome=chr,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=fit@parameters[n.segments, 1:n.segments])
  }
}
segs <- do.call(rbind, segs.list)
selection <- do.call(rbind, selection.list)
changes <- segs[1 < start,]
error.list <- labelError(
  selection, ann, changes,
  problem.vars="chromosome", # for all three data sets.
  model.vars="n.segments", # for changes and selection.
  change.var="chromStart", # column of changes with breakpoint position.
  label.vars=c("min", "max")) # limit of labels in ann.
pro.with.ann <- data.table(pro)[chromosome \%in\% ann$chromosome, ]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(n.segments ~ chromosome, scales="free", space="free")+
  scale_x_continuous(breaks=c(100, 200))+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                                 "false negative"=3,
                                 "false positive"=1))+
  scale_fill_manual("label", values=change.colors)+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6),
                color="grey",
                fill=NA,
                data=error.list$label.errors)+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6,
                    fill=annotation, linetype=status),
                data=error.list$label.errors)+
  geom_point(aes(position/1e6, logratio),
             data=pro.with.ann,
             shape=1)+
  geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
               data=segs,
               color="green",
               size=1)+
  geom_vline(aes(xintercept=chromStart/1e6),
             data=changes,
             linetype="dashed",
             color="green")
## The BIC model selection criterion is lambda = log(n), where n is
## the number of data points to segment. This implies log(lambda) =
## log(log(n)) = the log2.n feature in all.features.mat.
pred <- pro.with.ann[, list(pred.log.lambda=log(log(.N))), by=chromosome]
result <- ROChange(error.list$model.errors, pred, "chromosome")

ggplot()+
  geom_path(aes(FPR, TPR), data=result$roc)+
  geom_point(aes(FPR, TPR, color=threshold), data=result$thresholds, shape=1)

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
