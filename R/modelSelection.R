### Describe an annotated region label for supervised change-point detection.
changeLabel <- function(annotation, min.changes, max.changes, color){
  data.table(annotation, min.changes, max.changes, color)
}
change.labels <- rbind(
  changeLabel("breakpoint", 1, Inf, "#a445ee"),
  changeLabel(">0breakpoints", 1, Inf, "#a445ee"),
  changeLabel(">0changes", 1, Inf, "#a445ee"),
  changeLabel("normal", 0, 0, "#f6f4bf"),
  changeLabel("0breakpoints", 0, 0, "#f6f4bf"),
  changeLabel("0changes", 0, 0, "#f6f4bf"),
  changeLabel("1breakpoint", 1, 1, "#ff7d7d"),
  changeLabel("1change", 1, 1, "#ff7d7d"))
with(change.labels, stopifnot(min.changes <= max.changes))
stopifnot(0 <= change.labels$min.changes)
change.labels$possible.fn <- ifelse(0 < change.labels$min.changes, 1, 0)
change.labels$possible.fp <- ifelse(Inf == change.labels$max.changes, 0, 1)
setkey(change.labels, annotation)
change.colors <- paste(change.labels$color)
names(change.colors) <- change.labels$annotation

modelSelectionC <- structure(function # Exact model selection function
### Given loss.vec L_i, model.complexity K_i, the model selection
### function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
### the solutions (i, min.lambda, max.lambda) with i being the
### solution for every lambda in (min.lambda, max.lambda). This
### function uses the linear time algorithm implemented in C code.
### This function is mostly meant for internal use -- it is instead
### recommended to use modelSelection.
 (loss.vec,
### numeric vector: loss L_i
  model.complexity,
### numeric vector: model complexity K_i
  model.id
### vector: indices i
  ){
  stopifnot(is.numeric(loss.vec))
  stopifnot(is.numeric(model.complexity))
  stopifnot(0 < diff(model.complexity))
  stopifnot(diff(loss.vec) < 0)
  stopifnot(length(loss.vec) == length(model.complexity))
  stopifnot(length(model.id) == length(model.complexity))
  n.models <- length(loss.vec)
  after.vec <- rep(-1L, n.models)
  lambda.vec <- rep(-1, n.models)
  result.list <- .C(
    "modelSelection_interface",
    loss.vec=as.double(loss.vec),
    model.complexity=as.double(model.complexity),
    n.models=as.integer(n.models),
    after.vec=as.integer(after.vec),
    lambda.vec=as.double(lambda.vec))
  is.out <- 0 < result.list$lambda.vec
  lambda.out <- result.list$lambda.vec[is.out]
  i <- c(n.models, result.list$after.vec[is.out]+1)
  min.lambda <- c(0, lambda.out)
  max.lambda <- c(lambda.out, Inf)
  data.frame(
    min.lambda,
    max.lambda,
    min.log.lambda = log(min.lambda),
    max.log.lambda = log(max.lambda),
    model.complexity = model.complexity[i],
    model.id=model.id[i],
    model.loss=loss.vec[i],
    row.names=model.id[i])
### data.frame with a row for each model that can be selected for at
### least one lambda value, and the following columns. (min.lambda,
### max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
### optimal penalty constants, on the original and log scale;
### model.complexity are the K_i values; model.id are the model
### identifiers (also used for row names); and model.loss are the C_i
### values.
},ex=function(){

  data(neuroblastoma, package="neuroblastoma", envir=environment())
  pro <- subset(neuroblastoma$profiles, profile.id==1 & chromosome=="X")
  max.segments <- 20
  fit <- cghseg:::segmeanCO(pro$logratio, Kmax=max.segments)
  seg.vec <- 1:max.segments
  exact.df <- modelSelectionC(fit$J.est, seg.vec, seg.vec)
  ## Solve the optimization using grid search.
  L.grid <- with(exact.df,{
    seq(min(max.log.lambda)-1,
        max(min.log.lambda)+1,
        l=100)
  })
  lambda.grid <- exp(L.grid)
  kstar.grid <- sapply(lambda.grid, function(lambda){
    crit <- with(exact.df, model.complexity * lambda + model.loss)
    picked <- which.min(crit)
    exact.df$model.id[picked]
  })
  grid.df <- data.frame(log.lambda=L.grid, segments=kstar.grid)
  library(ggplot2)
  ## Compare the results.
  ggplot()+
    ggtitle("grid search (red) agrees with exact path computation (black)")+
    geom_segment(aes(min.log.lambda, model.id,
                     xend=max.log.lambda, yend=model.id),
                 data=exact.df)+
    geom_point(aes(log.lambda, segments),
               data=grid.df, color="red", pch=1)+
    ylab("optimal model complexity (segments)")+
    xlab("log(lambda)")
  
})

modelSelectionR <- structure(function # Exact model selection function
### Given loss.vec L_i, model.complexity K_i, the model selection
### function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
### the solutions (i, min.lambda, max.lambda) with i being the
### solution for every lambda in (min.lambda, max.lambda). This
### function uses the quadratic time algorithm implemented in R code.
### This function is mostly meant for internal use -- it is instead
### recommended to use modelSelection.
 (loss.vec,
### numeric vector: loss L_i
  model.complexity,
### numeric vector: model complexity K_i
  model.id
### vector: indices i
  ){
  stopifnot(is.numeric(loss.vec))
  stopifnot(is.numeric(model.complexity))
  stopifnot(0 < diff(model.complexity))
  stopifnot(diff(loss.vec) < 0)
  stopifnot(length(loss.vec) == length(model.complexity))
  stopifnot(length(model.id) == length(model.complexity))
  n.models <- length(loss.vec)
  Kmax <- model.complexity[n.models]
  Kcurrent <- Kmax
  Lcurrent <- 0
  vK <- Kmax
  vL <- 0
  vP <- model.id[n.models]
  i <- 2
  min.complexity <- model.complexity[1]
  while(Kcurrent > min.complexity) {
    is.smaller <- model.complexity < Kcurrent
    is.current <- model.complexity == Kcurrent
    smallerK <- model.complexity[is.smaller]
    smallerID <- model.id[is.smaller]
    loss.term <- loss.vec[is.current] - loss.vec[is.smaller]
    complexity.term <- smallerK - model.complexity[is.current]
    lambdaTransition <- loss.term/complexity.term
    next.i <- which.min(lambdaTransition)
    Kcurrent <- smallerK[next.i]
    Lcurrent <- min(lambdaTransition)
    vL[i] <- Lcurrent
    vK[i] <- Kcurrent
    vP[i] <- smallerID[next.i]
    i <- i + 1
  }
  L <- log(vL)
  data.frame(
    min.lambda = vL,
    max.lambda = c(vL[-1], Inf),
    min.log.lambda = L,
    max.log.lambda = c(L[-1], Inf),
    model.complexity = vK,
    model.id=vP,
    model.loss=loss.vec[vP],
    row.names=vP)
### data.frame with a row for each model that can be selected for at
### least one lambda value, and the following columns. (min.lambda,
### max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
### optimal penalty constants, on the original and log scale;
### model.complexity are the K_i values; model.id are the model
### identifiers (also used for row names); and model.loss are the C_i
### values.
},ex=function(){

  library(penaltyLearning)
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  one <- subset(neuroblastoma$profiles, profile.id==599 & chromosome=="14")
  max.segments <- 1000
  fit <- Segmentor3IsBack::Segmentor(one$logratio, model=2, Kmax=max.segments)
  lik.df <- data.frame(lik=fit@likelihood, segments=1:max.segments)
  times.list <- list()
  for(n.segments in seq(10, max.segments, by=10)){
    some.lik <- lik.df[1:n.segments,]
    some.times <- microbenchmark::microbenchmark(
      R=pathR <- with(some.lik, modelSelectionR(lik, segments, segments)),
      C=pathC <- with(some.lik, modelSelectionC(lik, segments, segments)),
      times=5)
    times.list[[paste(n.segments)]] <- data.frame(n.segments, some.times)
  }
  times <- do.call(rbind, times.list)
  ## modelSelectionR and modelSelectionC should give identical results.
  identical(pathR, pathC)
  ## However, modelSelectionC is much faster (linear time complexity)
  ## than modelSelectionR (quadratic time complexity).
  library(ggplot2)
  ggplot()+
    geom_point(aes(n.segments, time/1e9, color=expr), data=times)

})

largestContinuousMinimumR <- structure(function
### Find the run of minimum cost with the largest size.
(cost,
 size
 ){
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
  ##print(runs)
  largest <- which.max(runs$size)
  list(start=starts[largest],end=ends[largest])
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
  one.problem.error[indices$start:indices$end,]
  
})

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
    L <- largestContinuousMinimumR(errors, max.log.lambda-min.log.lambda)
    data.table(
      min.log.lambda=min.log.lambda[L$start],
      max.log.lambda=max.log.lambda[L$end])
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

labelError <- structure(function # Compute incorrect labels
### Compute incorrect labels for several change-point detection
### problems and models.
(models,
### data.frame with one row per (problem,model) combination.
  labels,
### data.frame with one row per (problem,label).
  changes,
### data.frame with one row per (problem,model,change).
  change.var="chromStart",
### character(length=1): column name of predicted change-point
### position (refers to the changes argument). The default
### "chromStart" is useful for genomic data with segment start/end
### positions stored in columns named chromStart/chromEnd.
  label.vars=c("min", "max"),
### character(length=2): column names of start and end positions of
### labeled regions, in same units as change-point positions (refers
### to the labels argument). The default is c("min", "max").
  model.vars="n.segments",
### character: column names used to identify model complexity. The
### default "n.segments" is for change-point models such as in the
### Segmentor3IsBack and cghseg packages.
  problem.vars=character(0)
### character: column names used to identify data set / segmentation
### problem. 
){
  stopifnot(is.character(problem.vars))
  stopifnot(is.character(model.vars))
  stopifnot(is.character(change.var))
  stopifnot(is.character(label.vars))
  stopifnot(length(change.var)==1)
  stopifnot(length(label.vars)==2)
  stopifnot(is.data.frame(models))
  stopifnot(is.data.frame(labels))
  stopifnot(is.data.frame(changes))
  if(length(problem.vars)==0){
    stop("Need at least one column name in problem.vars")
  }
  if(length(model.vars)==0){
    stop("Need at least one column name in model.vars")
  }
  stopifnot(label.vars %in% names(labels))
  stopifnot(change.var %in% names(changes))
  stopifnot(problem.vars %in% names(changes))
  stopifnot(problem.vars %in% names(labels))
  stopifnot(problem.vars %in% names(models))
  stopifnot(model.vars %in% names(models))
  stopifnot(model.vars %in% names(changes))
  new.key <- paste0(change.var, ".after")
  if(new.key %in% names(changes)){
    stop("changes should not have a column named ", new.key)
  }
  labels.dt <- data.table(labels)
  setkey(labels.dt, annotation)
  labels.info <- change.labels[labels.dt]
  stopifnot(nrow(labels.info)==nrow(labels.dt))
  setkeyv(labels.info, problem.vars)
  models.dt <- data.table(models)
  setkeyv(models.dt, problem.vars)
  model.labels <- labels.info[models.dt, allow.cartesian=TRUE]
  changes.dt <- data.table(changes)
  changes.dt[[new.key]] <- changes.dt[[change.var]]+1
  changes.key <- c(problem.vars, model.vars, change.var, new.key)
  setkeyv(changes.dt, changes.key)
  labels.key <- c(problem.vars, model.vars, label.vars)
  setkeyv(model.labels, labels.key)
  over.dt <- foverlaps(changes.dt, model.labels, nomatch=0L)
  long.key <- c(
    labels.key, "annotation",
    "min.changes", "max.changes", "possible.fp", "possible.fn")
  setkeyv(over.dt, long.key)
  setkeyv(model.labels, long.key)
  changes.per.label <- over.dt[model.labels, list(
    pred.changes=.N
  ), by=.EACHI]
  changes.per.label[, fp := ifelse(max.changes < pred.changes, 1, 0)]
  changes.per.label[, fn := ifelse(pred.changes < min.changes, 1, 0)]
  changes.per.label[, status := ifelse(
    fp, "false positive", ifelse(
      fn, "false negative", "correct"))]
  setkeyv(models.dt, c(problem.vars, model.vars))
  setkeyv(changes.per.label, c(problem.vars, model.vars))
  error.totals <- changes.per.label[models.dt, list(
    possible.fp=sum(possible.fp),
    fp=sum(fp),
    possible.fn=sum(possible.fn),
    fn=sum(fn),
    labels=.N,
    errors=sum(fp+fn)),
    by=.EACHI][models.dt]
  list(model.errors=error.totals, label.errors=changes.per.label)
### list of two data.tables: label.errors has one row for every
### combination of models and labels, with status column that
### indicates whether or not that model commits an error in that
### particular label; model.errors has one row per row of models, with
### columns for computing error and ROC curves.
}, ex=function(){
  
  library(penaltyLearning)
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  pro4 <- subset(neuroblastoma$profiles, profile.id==4)
  ann4 <- subset(neuroblastoma$annotations, profile.id==4)
  label <- function(annotation, min, max){
    data.table(profile.id=4, chromosome="14", min, max, annotation)
  }
  ann <- rbind(
      ann4,
      label("1change", 70e6, 80e6),
      label("0changes", 20e6, 60e6))
  max.segments <- 5
  segs.list <- list()
  models.list <- list()
  for(chr in unique(ann$chromosome)){
    pro <- subset(pro4, chromosome==chr)
    fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
    model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
    models.list[[chr]] <- data.table(chromosome=chr, model.df)
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
  models <- do.call(rbind, models.list)
  
  changes <- segs[1 < start,]
  error.list <- labelError(
      models, ann, changes,
      problem.vars="chromosome", # for all three data sets.
      model.vars="n.segments", # for changes and selection.
      change.var="chromStart", # column of changes with breakpoint position.
      label.vars=c("min", "max")) # limit of labels in ann.
  
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
               data=subset(pro4, chromosome %in% ann$chromosome),
               shape=1)+
    geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
                 data=segs,
                 color="green",
                 size=1)+
    geom_vline(aes(xintercept=chromStart/1e6),
               data=changes,
               linetype="dashed",
               color="green")
  
})

modelSelection <- function
### Given loss.vec L_i, model.complexity K_i, the model selection
### function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
### the solutions (i, min.lambda, max.lambda) with i being the
### solution for every lambda in (min.lambda, max.lambda). This
### function uses the quadratic time algorithm implemented in R code.
(models,
### data.frame with one row per model. There must be at
### least two columns [[loss]] and [[complexity]], but there can
### also be other meta-data columns.
 loss="loss",
### character: column name of models to interpret as loss L_i.
 complexity="complexity"
### character: column name of models to interpret as complexity K_i.
){
  stopifnot(is.data.frame(models))
  stopifnot(1 < nrow(models))
  for(x in list(loss, complexity)){
    stopifnot(is.character(x))
    stopifnot(length(x)==1)
    stopifnot(x %in% names(models))
  }
  ord <- order(models[[complexity]], models[[loss]])
  sorted <- models[ord,]
  keep <- c(TRUE, diff(sorted$error) < 0)
  filtered <- sorted[keep, ]
  loss.vec <- filtered[[loss]]
  complexity.vec <- filtered[[complexity]]
  id.vec <- (1:nrow(models))[ord][keep]
  result <- modelSelectionC(loss.vec, complexity.vec, id.vec)
  is.new <- names(result) %in% c(
    "model.complexity", "model.loss", "model.id")
  data.frame(
    result[!is.new],
    models[result$model.id,],
    row.names=rownames(models)[result$model.id])
### data.frame with a row for each model that can be selected for at
### least one lambda value, and the following columns. (min.lambda,
### max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
### optimal penalty constants, on the original and log scale;
### the other columns (and rownames) are taken from models.
}
