ROChange <- structure(function # ROC curve for changepoints
### Compute a ROC curve for a penalty function.
(models,
### data.frame describing the number of incorrect labels as a function
### of log(lambda), with columns min.log.lambda, max.log.lambda, fp,
### fn, possible.fp, possible.fn, etc. This can be computed via
### labelError(modelSelection(...), ...)$model.errors -- see examples.
  predictions,
### data.frame with the predicted log(lambda) value for each
### segmentation problem.
  problem.vars=character()
### character: column names used to identify data set / segmentation
### problem. 
){
  pred <- data.table(predictions)
  err <- data.table(models)
  setkey(err, min.log.lambda)
  thresh.dt <- err[, {
    fp.diff <- diff(fp)
    fp.change <- fp.diff != 0
    fn.diff <- diff(fn)
    fn.change <- fn.diff != 0
    fp.dt <- if(any(fp.change))data.table(
      log.lambda=max.log.lambda[c(fp.change, FALSE)],
      fp=fp.diff[fp.change],
      fn=0)
    fn.dt <- if(any(fn.change))data.table(
      log.lambda=max.log.lambda[c(fn.change, FALSE)],
      fp=0,
      fn=fn.diff[fn.change])
    rbind(fp.dt, fn.dt)
  }, by=problem.vars]
  setkey(thresh.dt, log.lambda)
  total.dt <- err[, .SD[1,], by=problem.vars][, list(
    labels=sum(labels),
    possible.fp=sum(possible.fp),
    possible.fn=sum(possible.fn))]
  setkeyv(thresh.dt, problem.vars)
  setkeyv(pred, problem.vars)
  pred.with.thresh <- pred[thresh.dt]
  pred.with.thresh[, thresh := log.lambda - pred.log.lambda]
  setkey(pred.with.thresh, thresh)
  interval.dt <- pred.with.thresh[, data.table(
    total.dt,
    min.thresh=c(-Inf, thresh),
    max.thresh=c(thresh, Inf),
    fp=total.dt$possible.fp + cumsum(c(0, fp)),
    fn=cumsum(c(0, fn)))]
  interval.dt[, errors := fp+fn]
  interval.dt[, FPR := fp/possible.fp]
  interval.dt[, tp := possible.fn - fn]
  interval.dt[, TPR := tp/possible.fn]
  interval.dt[, error.percent := 100*errors/labels]
  interval.dt
}, ex=function(){

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
  pro.with.ann <- data.table(pro)[chromosome %in% ann$chromosome, ]
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
  roc <- ROChange(error.list$model.errors, pred, "chromosome")

  pred.thresh <- roc[min.thresh < 0 & 0 < max.thresh,]
  ggplot()+
    geom_path(aes(FPR, TPR), data=roc)+
    geom_point(aes(FPR, TPR), data=pred.thresh, shape=1)

  ggplot()+
    geom_segment(aes(
      min.thresh, errors,
      xend=max.thresh, yend=errors),
      data=roc)+
    geom_point(aes(0, errors), data=pred.thresh, shape=1)+
    xlab("log(penalty) constant added to BIC penalty")

})

