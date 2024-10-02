\name{ROChange}
\alias{ROChange}
\title{ROC curve for changepoints}
\description{Compute a Receiver Operating Characteristic curve for a penalty
function.}
\usage{ROChange(models, predictions, 
    problem.vars = character())}
\arguments{
  \item{models}{data.frame describing the number of incorrect labels as a function
of log(lambda), with columns min.log.lambda, max.log.lambda, fp,
fn, possible.fp, possible.fn, etc. This can be computed via
labelError(modelSelection(...), ...)$model.errors -- see examples.}
  \item{predictions}{data.frame with a column named pred.log.lambda, the predicted
log(penalty) value for each segmentation problem.}
  \item{problem.vars}{character: column names used to identify data set / segmentation
problem.}
}

\value{named list of results:
\item{roc}{a data.table with one row for each point on
the ROC curve}
\item{thresholds}{two rows of roc which correspond to the
predicted and minimal error thresholds}
\item{auc.polygon}{a data.table with one row for
each vertex of the polygon used to
compute AUC}
\item{auc}{numeric Area Under the ROC curve}
\item{aum}{numeric Area Under Min(FP,FN)}
\item{aum.grad}{data.table with one row for each prediction, and
columns hi/lo bound for the aum
generalized gradient.}}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

  library(penaltyLearning)
  library(data.table)

  data(neuroblastomaProcessed, envir=environment())
  ## Get incorrect labels data for one profile.
  pid <- 11
  pro.errors <- neuroblastomaProcessed$errors[
    profile.id==pid][order(chromosome, min.log.lambda)]
  dcast(pro.errors, n.segments ~ chromosome, value.var="errors")
  ## Get the feature that corresponds to the BIC penalty = log(n),
  ## meaning log(penalty) = log(log(n)).
  chr.vec <- paste(c(1:4, 11, 17))
  pid.names <- paste0(pid, ".", chr.vec)
  BIC.feature <- neuroblastomaProcessed$feature.mat[pid.names, "log2.n"]
  pred <- data.table(pred.log.lambda=BIC.feature, chromosome=chr.vec)
  ## edit one prediction so that it ends up having the same threshold
  ## as another one, to illustrate an aum sub-differential with
  ## un-equal lo/hi bounds.
  err.changes <- pro.errors[, {
    .SD[c(NA, diff(errors) != 0), .(min.log.lambda)]
  }, by=chromosome]
  (ch.vec <- err.changes[, structure(min.log.lambda, names=chromosome)])
  other <- "11"
  (diff.other <- ch.vec[[other]]-pred[other, pred.log.lambda, on=.(chromosome)])
  pred["1", pred.log.lambda := ch.vec[["1"]]-diff.other, on=.(chromosome)]
  pred["4", pred.log.lambda := 2, on=.(chromosome)]
  ch.vec[["1"]]-pred["1", pred.log.lambda, on=.(chromosome)]
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

  ## Plot area under Min(FP,FN).
  err.colors <- c(
    "fp"="red",
    "fn"="deepskyblue",
    "min.fp.fn"="black")
  err.sizes <- c(
    "fp"=3,
    "fn"=2,
    "min.fp.fn"=1)
  roc.tall <- melt(result$roc, measure.vars=names(err.colors))
  area.rects <- data.table(
    chromosome="total",
    result$roc[0<min.fp.fn])
  (gg.total <- ggplot()+
     geom_vline(
       xintercept=0,
       color="grey")+
     geom_rect(aes(
       xmin=min.thresh, xmax=max.thresh,
       ymin=0, ymax=min.fp.fn),
       data=area.rects,
       alpha=0.5)+
     geom_text(aes(
       min.thresh, min.fp.fn/2,
       label=sprintf(
         "Area Under Min(FP,FN)=\%.3f ",
         result$aum)),
       data=area.rects[1],
       hjust=1,
       color="grey50")+
     geom_segment(aes(
       min.thresh, value,
       xend=max.thresh, yend=value,
       color=variable, size=variable),
       data=data.table(chromosome="total", roc.tall))+
     scale_size_manual(values=err.sizes)+
     scale_color_manual(values=err.colors)+
     theme_bw()+
     theme(panel.grid.minor=element_blank())+
     scale_x_continuous(
       "Prediction threshold")+
     scale_y_continuous(
       "Incorrectly predicted labels",
       breaks=0:10))

  ## Add individual error curves.
  tall.errors <- melt(
    pro.errors[pred, on=.(chromosome)],
    measure.vars=c("fp", "fn"))
  gg.total+
    geom_segment(aes(
      min.log.lambda-pred.log.lambda, value,
      xend=max.log.lambda-pred.log.lambda, yend=value,
      size=variable, color=variable),
      data=tall.errors)+
    facet_grid(chromosome ~ ., scales="free", space="free")+
    theme(panel.spacing=grid::unit(0, "lines"))+
    geom_blank(aes(
      0, errors),
      data=data.table(errors=c(1.5, -0.5)))

  print(result$aum.grad)
  if(interactive()){#this can be too long for CRAN.
    ## Plot how Area Under Min(FP,FN) changes with each predicted value.
    aum.dt <- pred[, {
      data.table(log.pen=seq(0, 4, by=0.5))[, {
        chr <- paste(chromosome)
        new.pred.dt <- data.table(pred)
        new.pred.dt[chr, pred.log.lambda := log.pen, on=.(chromosome)]
        with(
          ROChange(pro.errors, new.pred.dt, "chromosome"),
          data.table(aum))
      }, by=log.pen]
    }, by=chromosome]
    bounds.dt <- melt(
      result$aum.grad,
      measure.vars=c("lo", "hi"),
      variable.name="bound",
      value.name="slope")[pred, on=.(chromosome)]
    bounds.dt[, intercept := result$aum-slope*pred.log.lambda]
    ggplot()+
      geom_abline(aes(
        slope=slope, intercept=intercept),
        size=1,
        data=bounds.dt)+
      geom_text(aes(
        2, 2, label=sprintf("directional derivatives = [\%d, \%d]", lo, hi)),
        data=result$aum.grad)+
      scale_color_manual(
        values=c(
          predicted="red",
          new="black"))+
      geom_point(aes(
        log.pen, aum, color=type),
        data=data.table(type="new", aum.dt))+
      geom_point(aes(
        pred.log.lambda, result$aum, color=type),
        shape=1,
        data=data.table(type="predicted", pred))+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_wrap("chromosome", labeller=label_both)+
      coord_equal()+
      xlab("New log(penalty) value for chromosome")+
      ylab("Area Under Min(FP,FN)
using new log(penalty) for this chromosome
and predicted log(penalty) for others")
  }

}
