### stop with an informative error if there are problems with the
### target matrix or predicted values.
check_target_pred <- function(target.mat, pred){
  if(!{
    is.matrix(target.mat) &&
    is.numeric(target.mat) &&
    all(!is.na(target.mat)) &&
    ncol(target.mat)==2 &&
    all(target.mat[,1] < target.mat[,2])
  }){
    stop("target.mat must be a numeric matrix with two columns and no missing entries (lower and upper limit of target interval)")
  }
  if(!{
    is.numeric(pred) &&
      all(is.finite(pred))
  }){
    stop("pred must be a numeric vector or matrix with neither missing nor infinite entries")
  }
  if(is.matrix(pred)){
    if(nrow(pred) != nrow(target.mat)){
      stop("nrow(pred) must be same as nrow(target.mat)")
    }
  }else{
    if(length(pred) != nrow(target.mat)){
      stop("length(pred) must be same as nrow(target.mat)")
    }
  }
  if(any(apply(!is.finite(target.mat), 1, all))){
    stop("each row of target.mat must have at least one finite limit")
  }
  nrow(target.mat)
### number of observations.
}

targetIntervalROC <- structure(function
### Compute a ROC curve using a target interval matrix. A prediction
### less than the lower limit is considered a false positive (penalty
### too small, too many changes), and a prediction greater than the
### upper limit is a false negative (penalty too large, too few
### changes). WARNING: this ROC curve is less detailed than the one
### you get from ROChange! Use ROChange if possible.
(target.mat,
### n x 2 numeric matrix: target intervals of log(penalty) values that
### yield minimal incorrect labels.
 pred
### numeric vector: predicted log(penalty) values.
){
  min.log.lambda <- max.log.lambda <- errors <- fp <- fn <-
    possible.fp <- possible.fn <- NULL
### The code above is to avoid CRAN NOTEs like
### targetIntervals: no visible binding for global variable
  n <- check_target_pred(target.mat, pred)
  if(length(pred) != nrow(target.mat)){
    stop("length(pred) must be same as nrow(target.mat)")
  }
  observation <- 1:n
  target.errors <- rbind(
    data.table(
      observation,
      min.log.lambda=-Inf,
      max.log.lambda=ifelse(
        is.finite(target.mat[,1]), target.mat[,1], target.mat[,2]),
      fp=ifelse(is.finite(target.mat[,1]), 1, 0),
      fn=0),
    data.table(
      observation,
      min.log.lambda=target.mat[,1],
      max.log.lambda=target.mat[,2],
      fp=0,
      fn=0)[is.finite(target.mat[,1]) & is.finite(target.mat[,2])],
    data.table(
      observation,
      min.log.lambda=ifelse(
        is.finite(target.mat[,2]), target.mat[,2], target.mat[,1]),
      max.log.lambda=Inf,
      fp=0,
      fn=ifelse(is.finite(target.mat[,2]), 1, 0)))
  target.errors[, errors := fp + fn]
  target.errors[, labels := 1]
  target.errors[, possible.fp := max(fp), by=observation]
  target.errors[, possible.fn := max(fn), by=observation]
  pred.dt <- data.table(observation, pred.log.lambda=as.numeric(pred))
  ROChange(target.errors, pred.dt, "observation")  
### list describing ROC curves, same as ROChange.
}, ex=function(){

  library(penaltyLearning)
  data(neuroblastomaProcessed, envir=environment())

  pid.vec <- c("1", "4")
  chr <- 2
  incorrect.labels <-
    neuroblastomaProcessed$errors[profile.id%in%pid.vec & chromosome==chr]
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
  
})

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
  check_target_pred(target.mat, pred)
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
  data(neuroblastomaProcessed, envir=environment())
  ## The BIC model selection criterion is lambda = log(n), where n is
  ## the number of data points to segment. This implies log(lambda) =
  ## log(log(n)), which is the log2.n feature.
  row.name.vec <- grep(
    "^(4|520)[.]",
    rownames(neuroblastomaProcessed$feature.mat),
    value=TRUE)
  feature.mat <- neuroblastomaProcessed$feature.mat[row.name.vec, ]
  target.mat <- neuroblastomaProcessed$target.mat[row.name.vec, ]
  pred.dt <- data.table(
    row.name=row.name.vec,
    target.mat,
    feature.mat[, "log2.n", drop=FALSE])
  pred.dt[, pred.log.lambda := log2.n ]
  pred.dt[, residual := targetIntervalResidual(
    cbind(min.L, max.L),
    pred.log.lambda)]
  library(ggplot2)
  limits.dt <- pred.dt[, data.table(
    log2.n,
    log.penalty=c(min.L, max.L),
    limit=rep(c("min", "max"), each=.N))][is.finite(log.penalty)]
  ggplot()+
    geom_abline(slope=1, intercept=0)+
    geom_point(aes(
      log2.n,
      log.penalty,
      fill=limit),
      data=limits.dt,
      shape=21)+
    geom_segment(aes(
      log2.n, pred.log.lambda,
      xend=log2.n, yend=pred.log.lambda-residual),
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
  min.log.lambda <- errors <- max.log.lambda <- NULL
### The code above is to avoid CRAN NOTEs like
### targetIntervals: no visible binding for global variable
  stopifnot(is.data.frame(models))
  stopifnot(is.character(problem.vars))
  stopifnot(problem.vars %in% names(models))
  if(!is.numeric(models[["errors"]])){
    stop("models$errors should be the number of incorrect labels")
  }
  error.dt <- data.table(models)
  setkey(error.dt, min.log.lambda)
  error.dt[, {
    L <- largestContinuousMinimumC(errors, max.log.lambda-min.log.lambda)
    data.table(
      min.log.lambda=min.log.lambda[L[["start"]]],
      max.log.lambda=max.log.lambda[L[["end"]]],
      errors=errors[L[["end"]]])
  }, by=problem.vars]
### data.table with columns problem.vars, one row for each
### segmentation problem. The "min.log.lambda", and "max.log.lambda"
### columns give the largest interval of log(penalty) values which
### results in the minimum incorrect labels for that problem. This can
### be used to create the target.mat parameter of the
### IntervalRegression* functions.
}, ex=function(){

  library(penaltyLearning)
  data(neuroblastomaProcessed, envir=environment())
  targets.dt <- targetIntervals(
    neuroblastomaProcessed$errors,
    problem.vars=c("profile.id", "chromosome"))
  
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

  library(penaltyLearning)
  data(neuroblastomaProcessed, envir=environment())
  one.problem.error <-
    neuroblastomaProcessed$errors[profile.id=="4" & chromosome=="1"]
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

  library(penaltyLearning)
  data(neuroblastomaProcessed, envir=environment())
  one.problem.error <-
    neuroblastomaProcessed$errors[profile.id=="4" & chromosome=="1"]
  indices <- one.problem.error[, largestContinuousMinimumC(
    errors, max.log.lambda-min.log.lambda)]
  one.problem.error[indices[["start"]]:indices[["end"]],]
  
})

