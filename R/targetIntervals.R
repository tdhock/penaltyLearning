targetIntervalROC <- function
### Compute a ROC curve using a target interval matrix. WARNING: this
### ROC curve is less detailed than the one you get from ROChange! Use
### ROChange if possible.
(target.mat,
### n x 2 numeric matrix: target intervals of log(penalty) values that
### yield minimal incorrect labels.
 pred
### numeric vector: predicted log(penalty) values.
 ){
  stopifnot(is.matrix(target.mat))
  stopifnot(is.numeric(target.mat))
  stopifnot(target.mat[,1] < target.mat[,2])
  stopifnot(is.numeric(pred))
  stop("Not implemented")
### list describing ROC curves, same as ROChange.
}

targetIntervalResidual <- structure(function
### Compute residual of predicted penalties with respect to target
### intervals. This function is useful for visualizing the errors in a
### plot of log(penalty) versus a feature.
(target.mat,
### n x 2 numeric matrix: target intervals of log(penalty) values that
### yield minimal incorrect labels.
 pred
### numeric vector: predicted log(penalty) values.
 ){
  stopifnot(is.matrix(target.mat))
  stopifnot(is.numeric(target.mat))
  stopifnot(target.mat[,1] < target.mat[,2])
  stopifnot(is.numeric(pred))
  pred.vec <- as.numeric(pred)
  ifelse(
    pred.vec < target.mat[, 1], pred.vec - target.mat[, 1], ifelse(
      target.mat[, 2] < pred.vec, pred.vec - target.mat[, 2], 0))
### numeric vector of n residuals. Predictions that are too high
### (above target.mat[,2]) get positive residuals (too few
### changepoints), and predictions that are too low (below
### target.mat[,1]) get negative residuals.
}, ex=function(){

  library(penaltyLearning)
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  pid <- c(4, 520)
  pro <- subset(neuroblastoma$profiles, profile.id %in% pid)
  ann <- subset(neuroblastoma$annotations, profile.id %in% pid)
  max.segments <- 20
  segs.list <- list()
  selection.list <- list()
  for(i in 1:nrow(ann)){
    a <- ann[i,]
    id.str <- paste(a$chromosome, a$profile.id)
    pro.chr <- subset(pro, chromosome==a$chromosome & profile.id==a$profile.id)
    fit <- Segmentor3IsBack::Segmentor(
      pro.chr$logratio, model=2, Kmax=max.segments)
    model.dt <- data.table(
      lik=as.numeric(fit@likelihood),
      n.segments=1:max.segments,
      loss=NA_real_)
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
      seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
      data.mean.vec <- rep(seg.mean.vec, end-start+1)
      residual.vec <- pro.chr$logratio-data.mean.vec
      model.dt[n.segments, loss := sum(residual.vec * residual.vec)]
      segs.list[[paste(id.str, n.segments)]] <- data.table(
        profile.id=a$profile.id,
        chromosome=a$chromosome,
        n.segments,
        start,
        end,
        chromStart,
        chromEnd,
        mean=seg.mean.vec)
    }
    selection.dt <- modelSelection(model.dt, complexity="n.segments")
    selection.list[[id.str]] <- data.table(
      profile.id=a$profile.id,
      chromosome=a$chromosome,
      selection.dt)
  }
  segs <- do.call(rbind, segs.list)
  selection <- do.call(rbind, selection.list)
  changes <- segs[1 < start,]
  error.list <- labelError(
    selection, ann, changes,
    problem.vars=c("profile.id", "chromosome"),
    model.vars="n.segments", 
    change.var="chromStart", 
    label.vars=c("min", "max")) 
  target.dt <- targetIntervals(
    error.list$model.errors, c("profile.id", "chromosome"))
  ## The BIC model selection criterion is lambda = log(n), where n is
  ## the number of data points to segment. This implies log(lambda) =
  ## log(log(n)), which is a feature we can compute:
  feature.dt <- data.table(pro)[, list(
    loglog.n=log(log(.N))
  ), by=list(profile.id, chromosome)]
  pred.dt <- feature.dt[target.dt, on=list(profile.id, chromosome)]
  pred.dt[, pred.log.lambda := loglog.n ]
  pred.dt[, residual := targetIntervalResidual(
    cbind(min.log.lambda, max.log.lambda),
    pred.log.lambda)]

  library(ggplot2)
  limits.dt <- pred.dt[, data.table(
    loglog.n,
    log.penalty=c(min.log.lambda, max.log.lambda),
    limit=rep(c("min", "max"), each=.N))][is.finite(log.penalty)]
  ggplot()+
    geom_abline(slope=1, intercept=0)+
    geom_point(aes(
      loglog.n,
      log.penalty,
      fill=limit),
      data=limits.dt,
      shape=21)+
    geom_segment(aes(
      loglog.n, pred.log.lambda,
      xend=loglog.n, yend=pred.log.lambda-residual),
      data=pred.dt,
      color="red")+
    scale_fill_manual(values=c(min="white", max="black"))
  
})

targetIntervals <- structure(function # Compute target intervals
### Compute target intervals of log(penalty) values that result in
### predicted changepoint models with minimum incorrect labels.
### Use this function after labelError, and before IntervalRegression*.
(models,
### data.table with columns errors, min.log.lambda, max.log.lambda,
### typically labelError()$model.errors.
  problem.vars
### character: column names used to identify data set / segmentation
### problem.
){
  stopifnot(is.data.frame(models))
  stopifnot(is.character(problem.vars))
  stopifnot(problem.vars %in% names(models))
  error.dt <- data.table(models)
  setkey(error.dt, min.log.lambda)
  error.dt[, {
    L <- largestContinuousMinimumC(errors, max.log.lambda-min.log.lambda)
    data.table(
      min.log.lambda=min.log.lambda[L[["start"]]],
      max.log.lambda=max.log.lambda[L[["end"]]])
  }, by=problem.vars]
### data.table with columns problem.vars, one row for each
### segmentation problem. The "min.log.lambda", and "max.log.lambda"
### columns give the largest interval of log(penalty) values which
### results in the minimum incorrect labels for that problem. This can
### be used to create the target.mat parameter of the
### IntervalRegression* functions.
}, ex=function(){
  
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  pro4 <- subset(neuroblastoma$profiles, profile.id==4)
  ann4 <- subset(neuroblastoma$annotations, profile.id==4)
  label <- function(annotation, min, max){
    data.frame(profile.id=4, chromosome="14", min, max, annotation)
  }
  ann <- rbind(
    ann4,
    label("1change", 70e6, 80e6),
    label("0changes", 20e6, 60e6))
  max.segments <- 20
  segs.list <- list()
  selection.list <- list()
  for(chr in unique(ann$chromosome)){
    pro <- subset(pro4, chromosome==chr)
    fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
    model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
    selection.df <- modelSelection(model.df, complexity="n.segments")
    selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
    for(n.segments in 1:max.segments){
      end <- fit@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro$position))
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
  targetIntervals(error.list$model.errors, "chromosome")
  
})

largestContinuousMinimumR <- structure(function
### Find the run of minimum cost with the largest size.
### This function uses a two pass R implementation,
### and is meant for internal use.
### Use targetIntervals for real data.
(cost,
### numeric vector of cost values.
  size
### numeric vector of interval size values.
){
  stopifnot(
    is.numeric(cost),
    is.numeric(size),
    length(cost)==length(size),
    0 < size)
  m <- min(cost)
  is.min <- cost == m
  d <- c(diff(c(FALSE,is.min,FALSE)))
  ##print(data.frame(cost=c(cost,NA),size=c(size,NA),diff=d))
  starts <- which(d==1)
  ends <- which(d==-1)-1
  runs <- data.frame(starts,ends)
  ##print(runs)
  runs$size <- sapply(seq_along(starts),function(i){
    sum(size[ starts[i]:ends[i] ])
  })
  if(1 < sum(runs$size==Inf)){
    c(start=1L, end=length(cost))
  }else{
    largest <- which.max(runs$size)
    c(start=starts[largest],end=ends[largest])
  }
### Integer vector length 2 (start and end of target interval relative
### to cost and size).
}, ex=function(){

  data(neuroblastoma, package="neuroblastoma", envir=environment())
  pro4 <- subset(neuroblastoma$profiles, profile.id==4)
  ann4 <- subset(neuroblastoma$annotations, profile.id==4)
  label <- function(annotation, min, max){
    data.frame(profile.id=4, chromosome="14", min, max, annotation)
  }
  ann <- rbind(
    ann4,
    label("1change", 70e6, 80e6),
    label("0changes", 20e6, 60e6))
  max.segments <- 20
  segs.list <- list()
  selection.list <- list()
  for(chr in unique(ann$chromosome)){
    pro <- subset(pro4, chromosome==chr)
    fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
    model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
    selection.df <- modelSelection(model.df, complexity="n.segments")
    selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
    for(n.segments in 1:max.segments){
      end <- fit@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro$position))
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

  one.problem.error <- error.list$model.errors[chromosome=="14", ]
  indices <- one.problem.error[, largestContinuousMinimumR(
    errors, max.log.lambda-min.log.lambda)]
  one.problem.error[indices[["start"]]:indices[["end"]],]
  
})

largestContinuousMinimumC <- structure(function
### Find the run of minimum cost with the largest size.
### This function use a linear time C implementation,
### and is meant for internal use.
### Use targetIntervals for real data.
(cost,
### numeric vector of cost values.
  size
### numeric vector of interval size values.
){
  stopifnot(
    is.numeric(cost),
    is.numeric(size),
    length(cost)==length(size),
    !is.na(cost),
    !is.na(size),
    0 < size)
  result <- .C(
    "largestContinuousMinimum_interface",
    n_data=length(cost),
    cost_vec=as.double(cost),
    size_vec=as.double(size),
    index_vec=as.integer(c(0,0)),
    NAOK=TRUE,
    PACKAGE="penaltyLearning")
  indices <- result$index_vec + 1L
  names(indices) <- c("start", "end")
  indices
### Integer vector length 2 (start and end of target interval relative
### to cost and size).
}, ex=function(){

  data(neuroblastoma, package="neuroblastoma", envir=environment())
  pro4 <- subset(neuroblastoma$profiles, profile.id==4)
  ann4 <- subset(neuroblastoma$annotations, profile.id==4)
  label <- function(annotation, min, max){
    data.frame(profile.id=4, chromosome="14", min, max, annotation)
  }
  ann <- rbind(
    ann4,
    label("1change", 70e6, 80e6),
    label("0changes", 20e6, 60e6))
  max.segments <- 20
  segs.list <- list()
  selection.list <- list()
  for(chr in unique(ann$chromosome)){
    pro <- subset(pro4, chromosome==chr)
    fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
    model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
    selection.df <- modelSelection(model.df, complexity="n.segments")
    selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
    for(n.segments in 1:max.segments){
      end <- fit@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
      (pro$position[data.before.change]+pro$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro$position))
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

  one.problem.error <- error.list$model.errors[chromosome=="14", ]
  indices <- one.problem.error[, largestContinuousMinimumC(
    errors, max.log.lambda-min.log.lambda)]
  one.problem.error[indices[["start"]]:indices[["end"]],]
  
})

