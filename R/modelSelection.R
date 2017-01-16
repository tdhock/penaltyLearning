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
    prefix.vars="chromosome", # for all three data sets.
    model.vars="n.segments", # for changes and selection.
    change.var="chromStart", # column of changes with breakpoint position.
    label.vars=c("min", "max")) # limit of labels in ann.
  
})

### WANT: output list of two data.tables: label.errors has one row for
### every combination of models and labels, with status column that
### indicates whether or not that model commits an error in that
### particular label; model.errors has one row per row of models, with
### columns for computing error and ROC curves.
labelError <- structure(function
(models,
  labels,
  changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  model.vars=character(0), 
  prefix.vars=character(0)
){
  stopifnot(is.character(prefix.vars))
  stopifnot(is.character(model.vars))
  stopifnot(is.character(change.var))
  stopifnot(is.character(label.vars))
  stopifnot(length(change.var)==1)
  stopifnot(length(label.vars)==2)
  stopifnot(is.data.frame(models))
  stopifnot(is.data.frame(labels))
  stopifnot(is.data.frame(changes))
  if(length(prefix.vars)==0){
    stop("Need at least one column name in prefix.vars")
  }
  if(length(model.vars)==0){
    stop("Need at least one column name in model.vars")
  }
  stopifnot(label.vars %in% names(labels))
  stopifnot(change.var %in% names(changes))
  stopifnot(prefix.vars %in% names(changes))
  stopifnot(prefix.vars %in% names(labels))
  stopifnot(prefix.vars %in% names(models))
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
  setkeyv(labels.info, prefix.vars)
  models.dt <- data.table(models)
  setkeyv(models.dt, prefix.vars)
  model.labels <- labels.info[models.dt, allow.cartesian=TRUE]
  changes.dt <- data.table(changes)
  changes.dt[[new.key]] <- changes.dt[[change.var]]+1
  changes.key <- c(prefix.vars, model.vars, change.var, new.key)
  setkeyv(changes.dt, changes.key)
  labels.key <- c(prefix.vars, model.vars, label.vars)
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
  setkeyv(models.dt, c(prefix.vars, model.vars))
  setkeyv(changes.per.label, c(prefix.vars, model.vars))
  error.totals <- changes.per.label[models.dt, list(
    possible.fp=sum(possible.fp),
    fp=sum(fp),
    possible.fn=sum(possible.fn),
    fn=sum(fn),
    labels=.N,
    errors=sum(fp+fn)),
    by=.EACHI][models.dt]
  list(model.errors=error.totals, label.errors=changes.per.label)
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
    prefix.vars="chromosome", # for all three data sets.
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
