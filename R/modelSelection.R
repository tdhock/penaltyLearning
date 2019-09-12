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
  result.list <- .C(
    "modelSelectionLinear_interface",
    loss.vec=as.double(loss.vec),
    model.complexity=as.double(model.complexity),
    n.models=as.integer(n.models),
    selected.model.vec=integer(n.models),
    pen.break.vec=double(n.models),
    loop.eval.vec=integer(n.models),
    PACKAGE="penaltyLearning")
  index.vec <- (result.list$n.models+1L):1
  selected <- result.list$selected.model.vec[index.vec]+1L
  max.lambda <- result.list$pen.break.vec[index.vec]
  min.lambda <- c(0, max.lambda[-length(max.lambda)])
  data.frame(
    min.lambda,
    max.lambda,
    min.log.lambda = log(min.lambda),
    max.log.lambda = log(max.lambda),
    cum.iterations=cumsum(result.list$loop.eval.vec)[selected],
    model.complexity = model.complexity[selected],
    model.id=model.id[selected],
    model.loss=loss.vec[selected],
    row.names=model.id[selected])
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
  pro <- subset(neuroblastoma$profiles, profile.id==1 & chromosome=="X")
  max.segments <- 20
  fit <- Segmentor3IsBack::Segmentor(pro$logratio, 2, max.segments)
  seg.vec <- 1:max.segments
  exact.df <- modelSelectionC(fit@likelihood, seg.vec, seg.vec)
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
  iterations.vec <- 0
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
    iterations.vec[i] <- length(lambdaTransition)
    vP[i] <- smallerID[next.i]
    i <- i + 1
  }
  L <- log(vL)
  data.frame(
    min.lambda = vL,
    max.lambda = c(vL[-1], Inf),
    min.log.lambda = L,
    max.log.lambda = c(L[-1], Inf),
    cum.iterations=cumsum(iterations.vec),
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

  if(interactive()){
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
    library(data.table)
    times.dt <- data.table(times)
    times.dt[, seconds := time/1e9]
    gg <- ggplot()+
      geom_point(aes(
        n.segments, seconds, color=expr),
        data=times.dt,
        shape=1)
    gg+
      scale_x_log10()+
      scale_y_log10()
    R.times <- times.dt[expr=="R"]
    R.times[, log.seconds := log(seconds)]
    R.times[, log.segs := log(n.segments)]
    R.times[, segs.squared := n.segments*n.segments]
    thresh <- 700
    big.segs <- R.times[thresh < n.segments]
    (log.fit <- lm(log.seconds ~ log.segs, data=big.segs))
    (quad.fit <- lm(seconds ~ n.segments + segs.squared, data=R.times))
    big.segs[, pred.seconds := exp(predict(log.fit))]
    R.times[, pred.seconds := predict(quad.fit)]
    pred.dt <- rbind(
      big.segs[, data.table(model="log", expr="R", pred.seconds, n.segments)],
      R.times[, data.table(model="quad", expr="R", pred.seconds, n.segments)])
    ggfit <- gg+
      geom_line(aes(
        n.segments, pred.seconds, linetype=model),
        data=pred.dt)
    ggfit+
      scale_x_log10()+
      scale_y_log10()
    extrapolate.dt <- data.table(
      n.segments=c(100000, 250000, 500000, 1e6))
    extrapolate.dt[, log.segs := log(n.segments)]
    extrapolate.dt[, segs.squared := n.segments*n.segments]
    pred.seconds.mat <- rbind(
      quad=predict(quad.fit, extrapolate.dt),
      log=exp(predict(log.fit, extrapolate.dt)))
    colnames(pred.seconds.mat) <- extrapolate.dt$n.segments
    (pred.hours.mat <- pred.seconds.mat/60/60)
  }

})

modelSelection <- function # Compute exact model selection function
### Given loss.vec L_i, model.complexity K_i, the model selection
### function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
### the solutions (i, min.lambda, max.lambda) with i being the
### solution for every lambda in (min.lambda, max.lambda). Use this
### function after having computed changepoints and loss values for
### each model, and before using labelError. This function uses the
### linear time algorithm implemented in C code (modelSelectionC).
(models,
### data.frame with one row per model. There must be at least two
### columns models[[loss]] and models[[complexity]], but there can
### also be other meta-data columns.
 loss="loss",
### character: column name of models to interpret as loss L_i.
 complexity="complexity"
### character: column name of models to interpret as complexity K_i.
){
  if(!(
    is.character(complexity) &&
    length(complexity)==1
  )){
    stop("complexity must be a column name of models")
  }
  if(!(
    is.character(loss) &&
    length(loss)==1
  )){
    stop("loss must be a column name of models")
  }
  if(!(
    is.data.frame(models) &&
    0 < nrow(models) &&
    is.numeric(models[[complexity]]) &&
    is.numeric(models[[loss]]) &&
    all(!is.na(models[[complexity]])) &&
    all(!is.na(models[[loss]]))
  )){
    stop("models must be data.frame with at least one row and numeric columns models[[complexity]] and models[[loss]] which are not missing/NA")
  }
  ord <- order(models[[complexity]], models[[loss]])
  sorted <- models[ord,]
  keep <- c(TRUE, diff(sorted[[loss]]) < 0)
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
### optimal penalty constants, on the original and log scale; the
### other columns (and rownames) are taken from models. This should be
### used as the models argument of labelError.
}
