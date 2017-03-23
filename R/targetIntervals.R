targetIntervalResidual <- function
### Compute residual of predicted penalties with respect to target
### intervals.
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
}

targetIntervals <- structure(function # Compute target intervals
### Compute target intervals of log(penalty) values that result in
### predicted changepoint models with minimum incorrect labels.
(models,
### data.table with columns errors, min.log.lambda, max.log.lambda
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

