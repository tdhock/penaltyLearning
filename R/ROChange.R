ROChange <- structure(function # ROC curve for changepoints
### Compute a Receiver Operating Characteristic curve for a penalty
### function.
(models,
### data.frame describing the number of incorrect labels as a function
### of log(lambda), with columns min.log.lambda, max.log.lambda, fp,
### fn, possible.fp, possible.fn, etc. This can be computed via
### labelError(modelSelection(...), ...)$model.errors -- see examples.
  predictions,
### data.frame with a column named pred.log.lambda, the predicted
### log(lambda) value for each segmentation problem.
  problem.vars=character()
### character: column names used to identify data set / segmentation
### problem. 
){
  possible.fp <- possible.fn <- min.log.lambda <- fp <- fn <- thresh <-
    log.lambda <- pred.log.lambda <- errors <- FPR <- tp <- TPR <-
      error.percent <- min.thresh <- max.thresh <- NULL
### The code above is to avoid CRAN NOTEs like
### ROChange: no visible binding for global variable
  if(!(
    is.character(problem.vars) &&
      0 < length(problem.vars) &&
      !is.na(problem.vars)
    )){
    stop("problem.vars should be a character vector of column names (IDs for predictions and models)")
  }
  exp.cols <- c(
    "fp", "possible.fp", "fn", "possible.fn", "errors", "labels",
    problem.vars, "min.log.lambda")
  if(!(
    is.data.frame(models) &&
      all(exp.cols %in% names(models))
    )){
    stop("models should have columns ", paste(exp.cols, collapse=", "))
  }
  pred <- data.table(predictions)
  err <- data.table(models)
  total.dt <- err[pred, on=problem.vars, .SD[1,], by=problem.vars][, list(
    labels=sum(labels),
    possible.fp=sum(possible.fp),
    possible.fn=sum(possible.fn))]
  if(is.na(total.dt$labels)){
    stop("some predictions do not exist in models")
  }
  if(total.dt$possible.fp==0){
    stop("no negative labels")
  }
  if(total.dt$possible.fn==0){
    stop("no positive labels")
  }
  thresh.dt <- err[order(-min.log.lambda), {
    fp.diff <- diff(fp)
    fp.change <- fp.diff != 0
    fn.diff <- diff(fn)
    fn.change <- fn.diff != 0
    fp.dt <- if(any(fp.change))data.table(
      log.lambda=min.log.lambda[c(fp.change, FALSE)],
      fp=as.numeric(fp.diff[fp.change]),
      fn=0)
    fn.dt <- if(any(fn.change))data.table(
      log.lambda=min.log.lambda[c(fn.change, FALSE)],
      fp=0,
      fn=as.numeric(fn.diff[fn.change]))
    ##browser(expr=sample.id=="McGill0322")
    rbind(fp.dt, fn.dt)
  }, by=problem.vars]
  pred.with.thresh <- thresh.dt[pred, on=problem.vars]
  pred.with.thresh[, thresh := log.lambda - pred.log.lambda]
  uniq.thresh <- pred.with.thresh[, list(
    fp=sum(fp),
    fn=sum(fn)
  ), by=thresh]
  interval.dt <- uniq.thresh[order(-thresh), data.table(
    total.dt,
    min.thresh=c(thresh, -Inf),
    max.thresh=c(Inf, thresh),
    fp=cumsum(c(0, fp)),
    fn=total.dt$possible.fn+cumsum(c(0, fn)))]
  interval.dt[, errors := fp+fn]
  interval.dt[, FPR := fp/possible.fp]
  interval.dt[, tp := possible.fn - fn]
  interval.dt[, TPR := tp/possible.fn]
  interval.dt[, error.percent := 100*errors/labels]
  roc.polygon <- interval.dt[, {
    has11 <- FPR[.N]==1 & TPR[.N]==1
    has00 <- FPR[1]==0 & TPR[1]==0
    list(
      FPR=c(if(!has00)0, FPR, if(!has11)1, 1),
      TPR=c(if(!has00)0, TPR, if(!has11)1, 0)
      )
  }]
  list(
    roc=interval.dt,
    thresholds=rbind(
      data.table(
        threshold="predicted",
        interval.dt[min.thresh <= 0 & 0 <= max.thresh, ]),
      data.table(threshold="min.error", interval.dt[which.min(errors), ])),
    auc.polygon=roc.polygon,
    auc=roc.polygon[, geometry::polyarea(FPR, TPR)]
    )
### list of results describing ROC curve: roc is a data.table with one
### row for each point on the ROC curve; thresholds is the two rows of
### roc which correspond to the predicted and minimal error
### thresholds; auc.polygon is a data.table with one row for each
### vertex of the polygon used to compute AUC; auc is the numeric Area
### Under the ROC curve, actually computed via geometry::polyarea as
### the area inside the auc.polygon.
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
  library(ggplot2)
  ggplot()+
    theme_bw()+
    theme_no_space+
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

})

