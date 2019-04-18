### The squared hinge loss.
squared.hinge <- function(x, e=1){
  ifelse(x<e,(x-e)^2,0)
}

IntervalRegressionCVmargin <- structure(function
### Use cross-validation to fit an L1-regularized linear interval
### regression model by optimizing both margin and regularization
### parameters. This function just calls IntervalRegressionCV with a
### margin.vec parameter that is computed based on the finite target
### interval limits. If default parameters are used, this function
### should be about 10 times slower than IntervalRegressionCV
### (since this function computes n.margin=10 models
### per regularization parameter whereas IntervalRegressionCV
### only computes one).
### On large (N > 1000 rows) data sets,
### this function should yield a model which is a little
### more accurate than IntervalRegressionCV
### (since the margin parameter is optimized).
(feature.mat,
### Numeric feature matrix, n observations x p features.
  target.mat,
### Numeric target matrix, n observations x 2 limits.
  log10.diff=2,
### Numeric scalar: factors of 10 below the largest finite limit
### difference to use as a minimum margin value (difference on the
### log10 scale which is used to generate margin parameters). Bigger
### values mean a grid of margin parameters with a larger range. For
### example if the largest finite limit in target.mat is 26 and the
### smallest finite limit is -4 then the largest limit difference is
### 30, which will be used as the maximum margin parameter. If
### log10.diff is the default of 2 then that means the smallest margin
### parameter will be 0.3 (two factors of 10 smaller than 30).
  n.margin=10L,
### Integer scalar: number of margin parameters, by default 10.
  ...
### Passed to IntervalRegressionCV.
) {
  if(!(
    is.numeric(log10.diff) &&
    length(log10.diff)==1 &&
    is.finite(log10.diff)
  )){
    stop(paste(
      "log10.diff must be a finite numeric value",
      "(the number of factors of 10 below",
      "the largest limit difference",
      "to use as a minimum margin value)"))
  }
  if(!(
    is.integer(n.margin) &&
    length(n.margin)==1 &&
    is.finite(n.margin) &&
    0 < n.margin
  )){
    stop(paste(
      "n.margin must be a positive integer",
      "(number of margin parameters)"))
  }
  n.observations <- check_features_targets(feature.mat, target.mat)
  t.vec <- sort(target.mat[is.finite(target.mat)])
  d.vec <- diff(t.vec)
  pos.vec <- d.vec[0 < d.vec]
  log10.range <- log10(t.vec[length(t.vec)]-t.vec[1])
  IntervalRegressionCV(
    feature.mat,
    target.mat,
    ...,
    margin.vec=10^seq(
      log10.range-log10.diff,
      log10.range,
      l=n.margin))
### Model fit list from IntervalRegressionCV.
}, ex=function() {
  if(interactive()){
    library(penaltyLearning)
    data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
    if(require(future)){
      plan(multiprocess)
    }
    set.seed(1)
    fit <- with(neuroblastomaProcessed, IntervalRegressionCVmargin(
      feature.mat, target.mat, verbose=1))
    plot(fit)
    print(fit$plot.heatmap)
  }
})

IntervalRegressionCV <- structure(function
### Use cross-validation to fit an L1-regularized linear interval
### regression model by optimizing margin and/or regularization
### parameters.
### This function repeatedly calls IntervalRegressionRegularized, and by
### default assumes that margin=1. To optimize the margin,
### specify the margin.vec parameter
### manually, or use IntervalRegressionCVmargin
### (which takes more computation time
### but yields more accurate models).
### If the future package is available,
### two levels of future_lapply are used
### to parallelize on validation.fold and margin.
(feature.mat,
### Numeric feature matrix, n observations x p features.
  target.mat,
### Numeric target matrix, n observations x 2 limits.
  n.folds=ifelse(nrow(feature.mat) < 10, 3L, 5L),
### Number of cross-validation folds.
  fold.vec=sample(rep(1:n.folds, l=nrow(feature.mat))),
### Integer vector of fold id numbers.
  verbose=0,
### numeric: 0 for silent, bigger numbers (1 or 2) for more output.
  min.observations=10,
### stop with an error if there are fewer than this many observations.
  reg.type="min",
### Either "1sd" or "min" which specifies how the regularization
### parameter is chosen during the internal cross-validation
### loop. min: first take the mean of the K-CV error functions, then
### minimize it (this is the default since it tends to yield the least
### test error). 1sd: take the most regularized model with the same
### margin which is within one standard deviation of that minimum
### (this model is typically a bit less accurate, but much less
### complex, so better if you want to interpret the coefficients).
  incorrect.labels.db=NULL,
### either NULL or a data.table, which specifies the error function to
### compute for selecting the regularization parameter on the
### validation set. NULL means to minimize the squared hinge loss,
### which measures how far the predicted log(penalty) values are from
### the target intervals. If a data.table is specified, its first key
### should correspond to the rownames of feature.mat, and columns
### min.log.lambda, max.log.lambda, fp, fn, possible.fp, possible.fn;
### these will be used with ROChange to compute the AUC for each
### regularization parameter, and the maximimum will be selected (in
### the plot this is negative.auc, which is minimized). This
### data.table can be computed via
### labelError(modelSelection(...),...)$model.errors -- see
### example(ROChange). In practice this makes the computation longer,
### and it should only result in more accurate models if there are
### many labels per data sequence.
  initial.regularization=0.001,
### Passed to IntervalRegressionRegularized.
  margin.vec=1,
### numeric vector of margin size hyper-parameters. The computation
### time is linear in the number of elements of margin.vec -- more
### values takes more computation time, but yields slightly more
### accurate models (if there is enough data).
 ...
### passed to IntervalRegressionRegularized.
){
  validation.fold <- negative.auc <- threshold <- incorrect.labels <-
    variable <- value <- regularization <- folds <- status <- type <-
      vjust <- upper.limit <- lower <- upper <- fold <- NULL
### The code above is to avoid CRAN NOTEs like
### IntervalRegressionCV: no visible binding for global variable
  n.observations <- check_features_targets(feature.mat, target.mat)
  stopifnot(is.integer(n.folds))
  stopifnot(is.integer(fold.vec))
  stopifnot(length(fold.vec) == n.observations)
  stopifnot(
    is.character(reg.type),
    length(reg.type)==1,
    reg.type %in% c("1sd", "min"))
  if(n.observations < min.observations){
    stop(
      n.observations,
      " in data set but minimum of ",
      min.observations,
      "; decrease min.observations or use a larger data set")
  }
  all.finite <- apply(is.finite(feature.mat), 2, all)
  if(sum(all.finite)==0){
    stop("after filtering NA features, none remain for training")
  }
  if(!(
    is.numeric(margin.vec) &&
    0 < length(margin.vec) &&
    all(is.finite(margin.vec))
  )){
    stop("margin.vec must be a numeric vector of finite margin size parameters")
  }
  fold.limits <- data.table(
    lower=is.finite(target.mat[,1]),
    upper=is.finite(target.mat[,2]),
    fold=fold.vec)[, list(
      lower=sum(lower),
      upper=sum(upper)
    ), by=list(fold)]
  if(fold.limits[, any(lower==0 | upper==0)]){
    print(fold.limits)
    stop("some folds have no upper/lower limits; each fold should have at least one upper and one lower limit")
  }
  validation.fold.vec <- unique(fold.vec)
  LAPPLY <- if(requireNamespace("future.apply")){
    future.apply::future_lapply
  }else{
    lapply
  }
  validation.data.list <- LAPPLY(validation.fold.vec, function(validation.fold){
    ##print(validation.fold)
    is.validation <- fold.vec == validation.fold
    is.train <- !is.validation
    train.features <- feature.mat[is.train, all.finite, drop=FALSE]
    train.targets <- target.mat[is.train, , drop=FALSE]
    dt.list <- LAPPLY(margin.vec, function(margin){
      if(1 <= verbose){
        cat(sprintf("margin=%f vfold=%d\n", margin, validation.fold))
      }
      fit <- IntervalRegressionRegularized(
        train.features, train.targets, verbose=verbose,
        margin=margin,
        initial.regularization=initial.regularization
      )#, ...)
      validation.features <- feature.mat[is.validation, , drop=FALSE]
      pred.log.lambda <- fit$predict(validation.features)
      validation.targets <- target.mat[is.validation, , drop=FALSE]
      too.small <- pred.log.lambda < validation.targets[, 1]
      too.big <- validation.targets[, 2] < pred.log.lambda
      is.error <- too.small | too.big
      getLoss <- function(m)colMeans({
        squared.hinge(pred.log.lambda-validation.targets[, 1], m)+
        squared.hinge(validation.targets[, 2]-pred.log.lambda, m)
      })
      dt <- data.table(
        margin,
        validation.fold,
        regularization=fit$regularization.vec,
        squared.hinge.loss=getLoss(margin),
        mean.squared.error=getLoss(0),
        incorrect.intervals=colSums(is.error))
      if(!is.null(incorrect.labels.db)){
        dt$negative.auc <- NA_real_
        dt$incorrect.labels <- NA_real_
        for(regularization.i in seq_along(fit$regularization.vec)){
          pred.dt <- data.table(
            pred.log.lambda=pred.log.lambda[, regularization.i])
          pred.dt[[key(incorrect.labels.db)]] <-
            rownames(feature.mat)[is.validation]
          roc <- ROChange(
            incorrect.labels.db, pred.dt, key(incorrect.labels.db))
          dt[regularization.i, negative.auc := -roc$auc]
          predicted.thresh <- roc$thresholds[threshold=="predicted", ]
          dt[regularization.i, incorrect.labels := predicted.thresh$errors]
        }
      }
      dt
    })
    do.call(rbind, dt.list)
  })
  validation.data <- do.call(rbind, validation.data.list)
  vtall <- melt(
    validation.data,
    id.vars=c("validation.fold", "regularization", "margin"))
  vstats <- vtall[, list(
    mean=mean(value),
    sd=sd(value),
    folds=.N
  ), by=list(margin, regularization, variable)][folds==max(folds)]
  vstats.wide <- dcast(
    vstats, margin + regularization ~ variable, value.var="mean")
  validation.metrics <- if(is.null(incorrect.labels.db)){
    if(length(margin.vec)==1){
      ## only one margin parameter, so it is fine to use the squared
      ## hinge loss to choose the best model.
      "squared.hinge.loss"
    }else{
      ## several margin parameters, so it does not really make sense
      ## to compare them using the squared hinge loss. Also it does
      ## not make sense to use the mean squared error, since that
      ## favors models with small margin sizes. So we use the number
      ## of incorrect targets, and maybe the squared hinge loss to
      ## break ties.
      c("incorrect.intervals", "squared.hinge.loss")
    }
  }else{
    ## When we have AUC, use it first and then use other metrics to
    ## break ties.
    c("negative.auc",
      "incorrect.labels",
      "squared.hinge.loss")
  }
  ord.arg.list <- lapply(validation.metrics, function(N)vstats.wide[[N]])
  ord.vec <- do.call(order, ord.arg.list)
  validation.ord <- vstats.wide[ord.vec]
  validation.best <- validation.ord[1]
  best.margin.stats <- vstats[margin==validation.best$margin]
  first.variable <- best.margin.stats[variable==validation.metrics[1], ]
  best.regularization <-
    first.variable[regularization==validation.best$regularization]
  best.regularization[, upper.limit := mean+sd]
  within.1sd <- first.variable[mean < best.regularization$upper.limit]
  least.complex <- within.1sd[which.max(regularization)]
  dot.dt <- rbind(
    if(nrow(least.complex))data.table(
      type="1sd", least.complex, upper.limit=NA_real_),
    data.table(type="min", best.regularization))
  best.margin.folds <- vtall[margin==validation.best$margin]
  color.scale <- scale_color_manual(
    values=c(
      "1sd"="blue",
      "min"="red"))
  selected <- dot.dt[type==reg.type]
  if(nrow(selected)==0){
    stop(
      "reg.type=", reg.type,
      " undefined; try another reg.type (",
      paste(dot.dt$type, collapse=","),
      ") or decrease initial.regularization")
  }
  gg.bands <- ggplot()+
    ggtitle(paste0(
      "Regularization parameter selection using ",
      length(validation.fold.vec),
      "-fold cross-validation, margin=",
      validation.best$margin
    ))+
    theme_bw()+
    guides(color="none")+
    theme_no_space()+
    facet_grid(variable ~ ., scales="free")+
    color.scale+
    geom_vline(aes(
      xintercept=-log10(regularization),
      color=type
    ), data=selected)+
    geom_ribbon(aes(
      -log10(regularization),
      ymin=mean-sd,
      ymax=mean+sd),
      fill="grey",
      alpha=0.5,
      data=best.margin.stats)+
      geom_line(aes(
        -log10(regularization),
        value,
        group=validation.fold),
        color="grey",
        data=best.margin.folds)+
        geom_line(aes(
          -log10(regularization),
          mean),
          size=1,
          data=best.margin.stats)+
          geom_text(aes(-log10(regularization), mean,
                        label=paste0(type, " "),
                        color=type),
                    vjust=1,
                    hjust=1,
                    data=dot.dt)+
              geom_point(aes(-log10(regularization), mean,
                             color=type),
                         data=dot.dt)+
              xlab("model complexity -log(regularization)")+
              ylab("")
  gg.heatmap <- ggplot()+
    geom_tile(aes(
      -log10(regularization),
      log10(margin),
      fill=log10(mean)
    ), data=vstats[variable==validation.metrics[1],])+
    scale_fill_gradient(low="white", high="red")+
    color.scale+
    geom_point(aes(
      -log10(regularization),
      log10(margin),
      color=type
    ), data=dot.dt)
  fit <- IntervalRegressionRegularized(
    feature.mat, target.mat,
    initial.regularization=selected$regularization,
    factor.regularization=NULL,
    margin=selected$margin,
    verbose=verbose,
    ...)
  fit$plot.selectRegularization <- fit$plot <- gg.bands
  fit$plot.selectRegularization.line <- best.margin.folds
  fit$plot.selectRegularization.ribbon <- best.margin.stats
  fit$plot.selectRegularization.point <- dot.dt
  fit$plot.selectRegularization.vlines <- selected
  fit$plot.heatmap <- gg.heatmap
  fit$plot.heatmap.tile <- vstats
  fit$validation.data <- validation.data
  fit
### List representing regularized linear model.
}, ex=function(){

  if(interactive()){
    library(penaltyLearning)
    data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
    if(require(future)){
      plan(multiprocess)
    }
    set.seed(1)
    i.train <- 1:100
    fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
      feature.mat[i.train,], target.mat[i.train,],
      verbose=0))
    ## When only features and target matrices are specified for
    ## training, the squared hinge loss is used as the metric to
    ## minimize on the validation set.
    plot(fit)
    ## Create an incorrect labels data.table (first key is same as
    ## rownames of feature.mat and target.mat).
    library(data.table)
    errors.per.model <- data.table(neuroblastomaProcessed$errors)
    errors.per.model[, pid.chr := paste0(profile.id, ".", chromosome)]
    setkey(errors.per.model, pid.chr)
    set.seed(1)
    fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
      feature.mat[i.train,], target.mat[i.train,],
      ## The incorrect.labels.db argument is optional, but can be used if
      ## you want to use AUC as the CV model selection criterion.
      incorrect.labels.db=errors.per.model))
    plot(fit)
  }

})

IntervalRegressionUnregularized <- function
### Use IntervalRegressionRegularized with initial.regularization=0
### and factor.regularization=NULL, meaning fit one un-regularized
### interval regression model.
(...
### passed to IntervalRegressionRegularized.
 ){
  IntervalRegressionRegularized(
    ...,
    initial.regularization=0,
    factor.regularization=NULL)
### List representing fit model, see
### help(IntervalRegressionRegularized) for details.
}

check_features_targets <- function
### stop with an informative error if there is a problem with the
### feature or target matrix.
(feature.mat,
### n x p numeric input feature matrix.
  target.mat
### n x 2 matrix of target interval limits.
  ){
  if(!(
    is.numeric(feature.mat) &&
    is.matrix(feature.mat) &&
    is.character(colnames(feature.mat))
  )){
    stop("feature.mat should be a numeric matrix with colnames (input features)")
  }
  if(!(
    is.numeric(target.mat) &&
    is.matrix(target.mat) &&
    ncol(target.mat) == 2
  )){
    stop("target.mat should be a numeric matrix with two columns (lower and upper limits of correct outputs)")
  }
  if(!(
    any(is.finite(target.mat[,1]))
  )){
    stop("target.mat has no lower limits, but should have at least one")
  }
  if(!(
    any(is.finite(target.mat[,2]))
  )){
    stop("target.mat has no upper limits, but should have at least one")
  }
  if(nrow(target.mat) != nrow(feature.mat)){
    stop("feature.mat and target.mat should have the same number of rows")
  }
  nrow(feature.mat)
### number of observations/rows.
}

IntervalRegressionRegularized <- structure(function
### Repeatedly use IntervalRegressionInternal to solve interval
### regression problems for a path of regularization parameters. This
### function does not perform automatic selection of the
### regularization parameter; instead, it returns regression models
### for a range of regularization parameters, and it is up to you to
### select which one to use. For automatic regularization parameter
### selection, use IntervalRegressionCV.
(feature.mat,
### Numeric feature matrix.
 target.mat,
### Numeric target matrix.
 initial.regularization=0.001,
### Initial regularization parameter.
 factor.regularization=1.2,
### Increase regularization by this factor after finding an optimal
### solution. Or NULL to compute just one model
### (initial.regularization).
 verbose=0,
### Print messages if >= 1.
 margin=1,
### Non-negative margin size parameter, default 1.
 ...
### Other parameters to pass to IntervalRegressionInternal.
){
  residual <- limit <- normalized.weight <- variable <- NULL
### The code above is to avoid CRAN NOTEs like
### IntervalRegressionRegularized: no visible binding for global variable
  check_features_targets(feature.mat, target.mat)
  stopifnot(is.numeric(initial.regularization))
  stopifnot(length(initial.regularization)==1)
  stopifnot(is.finite(initial.regularization))
  is.trivial.target <- apply(!is.finite(target.mat), 1, all)
  nontrivial.features <- feature.mat[!is.trivial.target, , drop=FALSE]
  nontrivial.targets <- target.mat[!is.trivial.target, , drop=FALSE]
  is.finite.feature <- apply(is.finite(nontrivial.features), 2, all)
  finite.features <- nontrivial.features[, is.finite.feature, drop=FALSE]
  all.mean.vec <- colMeans(finite.features)
  all.sd.vec <- if(nrow(finite.features)==1){
    1
  }else{
    apply(finite.features, 2, sd)
  }
  is.invariant <- all.sd.vec == 0
  train.feature.i <- which(!is.invariant)
  train.feature.names <- colnames(finite.features)[train.feature.i]
  if(length(train.feature.names)==0){
    stop("after filtering NA and constant features, none remain for training")
  }
  mean.vec <- all.mean.vec[train.feature.names]
  sd.vec <- all.sd.vec[train.feature.names]
  invariant.features <- finite.features[, train.feature.names, drop=FALSE]
  mean.mat <- matrix(
    mean.vec, nrow(invariant.features), ncol(invariant.features), byrow=TRUE)
  sd.mat <- matrix(
    sd.vec, nrow(invariant.features), ncol(invariant.features), byrow=TRUE)
  norm.features <- (invariant.features-mean.mat)/sd.mat
  intercept.features <- cbind("(Intercept)"=1, norm.features)
  apply(intercept.features, 2, mean)
  apply(intercept.features, 2, sd)
  regularization <- initial.regularization
  n.features <- ncol(intercept.features)
  param.vec <- rep(0, n.features)
  n.nonzero <- n.features
  param.vec.list <- list()
  scaled.vec.list <- list()
  regularization.vec.list <- list()
  while(n.nonzero > 1){
    param.vec <-
      IntervalRegressionInternal(
        intercept.features, nontrivial.targets,
        param.vec,
        regularization,
        verbose=verbose,
        margin=margin,
        ...)
    n.zero <- sum(param.vec == 0)
    n.nonzero <- sum(param.vec != 0)
    l1.norm <- sum(abs(param.vec[-1]))
    if(verbose >= 1){
      cat(sprintf("regularization=%8.4f L1norm=%8.4f zeros=%d\n",
                  regularization, l1.norm, n.zero))
    }
    weight.vec <- param.vec[-1]
    ## training is done in the centered and scaled space
    ## (intercept.features), but we report coefficients for the
    ## original space.
    orig.param.vec <- c(
      param.vec[1] - sum(weight.vec*mean.vec/sd.vec),
      weight.vec/sd.vec)
    pred.vec <- intercept.features %*% param.vec
    orig.pred.vec <- cbind(1, invariant.features) %*% orig.param.vec
    stopifnot(all.equal(pred.vec, orig.pred.vec))
    scaled.vec.list[[paste(regularization)]] <- weight.vec
    param.vec.list[[paste(regularization)]] <- orig.param.vec
    regularization.vec.list[[paste(regularization)]] <- regularization
    if(is.null(factor.regularization)){
      n.nonzero <- 1 #stops while loop.
    }else{
      stopifnot(is.numeric(factor.regularization))
      stopifnot(length(factor.regularization)==1)
      regularization <- regularization * factor.regularization
    }
  }
  scaled.mat <- do.call(cbind, scaled.vec.list)
  param.mat <- do.call(cbind, param.vec.list)
  if(verbose >= 1){
    cat(paste0("Done computing parameter matrix (",
               nrow(param.mat), " features x ",
               ncol(param.mat), " regularization parameters)\n"))
  }
  feature.not.used <- apply(param.mat[-1, , drop=FALSE] == 0, 1, all)
  pred.feature.names <- train.feature.names[!feature.not.used]
  pred.param.mat <-
    param.mat[c("(Intercept)", pred.feature.names),,drop=FALSE]
  L <- list(
    margin=margin,
    param.mat=param.mat,
    regularization.vec=do.call(c, regularization.vec.list),
    train.feature.names=train.feature.names,
    pred.feature.names=pred.feature.names,
    pred.param.mat=pred.param.mat,
    plot.weight.data=data.table(
      normalized.weight=as.numeric(scaled.mat),
      original.weight=as.numeric(param.mat[-1,]),
      variable=rownames(scaled.mat)[row(scaled.mat)],
      regularization=as.numeric(colnames(scaled.mat)[col(scaled.mat)])),
    predict=function(mat){
      if(missing(mat))mat <- feature.mat
      stopifnot(is.matrix(mat))
      stopifnot(is.numeric(mat))
      is.missing <- ! pred.feature.names %in% colnames(mat)
      if(any(is.missing)){
        stop("columns needed for prediction but not present: ",
             paste(pred.feature.names[is.missing], collapse=", "))
      }
      cbind(1, mat[, pred.feature.names, drop=FALSE]) %*% pred.param.mat
    })
  class(L) <- c("IntervalRegression", "list")
  pred.log.penalty <- predict(L)
  lower.limit <- -Inf < target.mat[,1]
  upper.limit <- target.mat[,2] < Inf
  pred.dt <- data.table(
    pred.log.penalty=as.numeric(pred.log.penalty),
    regularization=as.numeric(
      colnames(pred.log.penalty)[col(pred.log.penalty)]),
    model.i=as.integer(col(pred.log.penalty)),
    observation=as.integer(row(pred.log.penalty)),
    residual=targetIntervalResidual(target.mat, pred.log.penalty),
    lower.limit, upper.limit, limit=ifelse(
      upper.limit, ifelse(lower.limit, "both", "upper"), "lower")
  )[lower.limit | upper.limit]
  L$plot.residual.data <- pred.dt
  limit.colors <- c(
    upper="red",
    both="blue",
    lower="grey50")
  L$plot.residual <- ggplot()+
    theme_bw()+
    theme_no_space()+
    facet_wrap("model.i")+
    geom_hline(yintercept=0, color="grey")+
    scale_color_manual(values=limit.colors, breaks=names(limit.colors))+
    geom_point(aes(
      pred.log.penalty, residual, color=limit),
      data=pred.dt,
      shape=1)
  L$plot <- if(1 == length(L$regularization.vec)){
    L$plot.residual
  }else{
    gg <- ggplot()+
      geom_line(aes(
        -log10(regularization), normalized.weight, color=variable),
        data=L$plot.weight.data)
    (L$plot.weight <- if(requireNamespace("directlabels")){
      directlabels::direct.label(gg, "lasso.labels")
    }else{
      message('install.packages("directlabels") for more informative labels on plot.weight')
      gg
    })
  }
  L
### List representing fit model. You can do
### fit$predict(feature.matrix) to get a matrix of predicted log
### penalty values. The param.mat is the n.features * n.regularization
### numeric matrix of optimal coefficients (on the original scale).
}, ex=function(){

  if(interactive()){
    library(penaltyLearning)
    data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
    i.train <- 1:500
    fit <- with(neuroblastomaProcessed, IntervalRegressionRegularized(
      feature.mat[i.train,], target.mat[i.train,]))
    plot(fit)
  }

})

### print learned model parameters.
print.IntervalRegression <- function(x, ...){
  if(ncol(x$pred.param.mat)==1){
    cat(
      "IntervalRegression model for margin=",
      x$margin, " regularization=",
      x$regularization.vec,
      " with weights:\n",
      sep="")
    x <- t(x$pred.param.mat)
    rownames(x) <- ""
    print(x)
  }else{
    cat(
      "IntervalRegression models for margin=",
      x$margin,
      " [",
      nrow(x$pred.param.mat),
      " weights x ",
      ncol(x$pred.param.mat),
      " regularizations]\n",
      sep="")
  }
}

### Get the learned coefficients of an IntervalRegression model.
coef.IntervalRegression <- function(object, ...){
  object$pred.param.mat
### numeric matrix [features x regularizations] of learned weights (on
### the original feature scale), can be used for prediction via
### cbind(1,features) %*% weights.
}

### Plot an IntervalRegression model.
plot.IntervalRegression <- function(x, ...){
  x$plot
### a ggplot.
}

### Compute model predictions.
predict.IntervalRegression <- function(object, X, ...){
  object$predict(X)
### numeric matrix of predicted log(penalty) values.
}

IntervalRegressionInternal <- function
### Solve the squared hinge loss interval regression problem for one
### regularization parameter: w* = argmin_w L(w) + regularization *
### ||w||_1 where L(w) is the average squared hinge loss with respect
### to the targets, and ||w||_1 is the L1-norm of the weight vector
### (excluding the first element, which is the un-regularized
### intercept or bias term). This function performs no scaling of
### input features, and is meant for internal use only! To learn a
### regression model, try IntervalRegressionCV or
### IntervalRegressionUnregularized.
(features,
### Scaled numeric feature matrix (problems x features). The first
### column/feature should be all ones and will not be regularized.
 targets,
### Numeric target matrix (problems x 2).
 initial.param.vec,
### initial guess for weight vector (features).
 regularization,
### Degree of L1-regularization.
 threshold=1e-3,
### When the stopping criterion gets below this threshold, the
### algorithm stops and declares the solution as optimal.
 max.iterations=1e3,
### If the algorithm has not found an optimal solution after this many
### iterations, increase Lipschitz constant and max.iterations.
 weight.vec=NULL,
### A numeric vector of weights for each training example.
 Lipschitz=NULL,
### A numeric scalar or NULL, which means to compute Lipschitz as the
### mean of the squared L2-norms of the rows of the feature matrix.
 verbose=2,
### Cat messages: for restarts and at the end if >= 1, and for every
### iteration if >= 2.
 margin=1
### Margin size hyper-parameter, default 1.
 ){
  if(!(
    is.numeric(margin) &&
    length(margin)==1 &&
    is.finite(margin)
    )){
    stop("margin must be finite numeric scalar")
  }
  stopifnot(is.matrix(features))
  stopifnot(is.numeric(features))
  n.features <- ncol(features)
  n.problems <- nrow(features)

  stopifnot(is.matrix(targets))
  stopifnot(nrow(targets) == n.problems)
  stopifnot(ncol(targets) == 2)

  if(is.null(weight.vec)){
    weight.vec <- rep(1, n.problems)
  }
  stopifnot(is.numeric(weight.vec))
  stopifnot(length(weight.vec) == n.problems)

  if(is.null(Lipschitz)){
    Lipschitz <- mean(rowSums(features * features) * weight.vec)
  }
  stopifnot(is.numeric(Lipschitz))
  stopifnot(length(Lipschitz) == 1)

  stopifnot(is.numeric(max.iterations))
  stopifnot(length(max.iterations) == 1)

  stopifnot(is.numeric(threshold))
  stopifnot(length(threshold) == 1)

  stopifnot(is.numeric(initial.param.vec))
  stopifnot(length(initial.param.vec) == n.features)

  ## Return 0 for a negative number and the same value otherwise.
  positive.part <- function(x){
    ifelse(x<0, 0, x)
  }
  squared.hinge.deriv <- function(x,e=1){
    ifelse(x<e,2*(x-e),0)
  }
  calc.loss <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge(linear.predictor-targets[,1], margin)
    right.term <- squared.hinge(targets[,2]-linear.predictor, margin)
    both.terms <- left.term+right.term
    weighted.loss.vec <- both.terms * weight.vec
    mean(weighted.loss.vec)
  }
  calc.grad <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge.deriv(linear.predictor-targets[,1], margin)
    right.term <- squared.hinge.deriv(targets[,2]-linear.predictor, margin)
    full.grad <- features * (left.term-right.term) * weight.vec
    colSums(full.grad)/nrow(full.grad)
  }
  calc.penalty <- function(x){
    regularization * sum(abs(x[-1]))
  }
  calc.cost <- function(x){
    calc.loss(x) + calc.penalty(x)
  }
  soft.threshold <- function(x,thresh){
    ifelse(abs(x) < thresh, 0, x-thresh*sign(x))
  }
  ## do not threshold the intercept.
  prox <- function(x,thresh){
    x[-1] <- soft.threshold(x[-1],thresh)
    x
  }
  ## p_L from the fista paper.
  pL <- function(x,L){
    grad <- calc.grad(x)
    prox(x - grad/L, regularization/L)
  }
  dist2subdiff.opt <- function(w,g){
    ifelse(w==0,positive.part(abs(g)-regularization),
           ifelse(w<0,abs(-regularization+g),abs(regularization+g)))
  }

  iterate.count <- 1
  stopping.crit <- threshold
  last.iterate <- this.iterate <- y <- initial.param.vec
  this.t <- 1
  while({
    ##browser(expr=is.na(stopping.crit))
    ##str(stopping.crit)
    stopping.crit >= threshold
  }){
    ## here we implement the FISTA method with constant step size, as
    ## described by in the Beck and Tebolle paper.
    last.iterate <- this.iterate
    this.iterate <- pL(y, Lipschitz)
    last.t <- this.t
    this.t <- (1+sqrt(1+4*last.t^2))/2
    y <- this.iterate + (last.t - 1)/this.t*(this.iterate-last.iterate)
    ## here we calculate the subgradient optimality condition, which
    ## requires 1 more gradient evaluation per iteration.
    after.grad <- calc.grad(this.iterate)
    w.dist <- dist2subdiff.opt(this.iterate[-1],after.grad[-1])
    zero.at.optimum <- c(abs(after.grad[1]),w.dist)
    stopping.crit <- max(zero.at.optimum)

    if(verbose >= 2){
      cost <- calc.cost(this.iterate)
      cat(sprintf("%10d cost %10f crit %10.7f\n",
                  iterate.count,
                  cost,
                  stopping.crit))
    }
    iterate.count <- iterate.count + 1
    if(any(!is.finite(this.iterate)) || 1e20 < stopping.crit){
      if(verbose >= 1){
        cat("restarting with bigger Lipschitz.\n")
      }
      iterate.count <- 1
      stopping.crit <- threshold
      last.iterate <- this.iterate <- y <- initial.param.vec
      this.t <- 1
      Lipschitz <- Lipschitz * 1.5
    }
    if(iterate.count > max.iterations){
      if(verbose >= 1){
        cat(max.iterations, "iterations, increasing Lipschitz and iterations.",
            "crit =", stopping.crit, "\n")
      }
      Lipschitz <- Lipschitz * 1.5
      iterate.count <- 1
      max.iterations <- max.iterations * 2
    }
  }
  if(verbose >= 1){
    cat("solution with crit =", stopping.crit, "\n")
  }
  this.iterate
### Numeric vector of scaled weights w of the affine function f_w(X) =
### X %*% w for a scaled feature matrix X with the first row entirely
### ones.
}

