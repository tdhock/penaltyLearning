### The squared hinge loss.
squared.hinge <- function(x){
  ifelse(x<1,(x-1)^2,0)
}

IntervalRegressionCV <- structure(function
### Use cross-validation to estimate the optimal regularization, by
### picking the value that minimizes the number of incorrectly
### predicted target intervals. K-fold cross-validation is
### parallelized using the foreach package.
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
 reg.type="min(mean)",
### Either "1sd", "min(mean)" or "mean(min)" which specifies how the
### regularization parameter is chosen during the internal
### cross-validation loop. min(mean): first take the mean of the K-CV
### error functions, then minimize it (this is the default since it
### tends to yield the least test error). 1sd: take the least complex
### model which is within one standard deviation of that minimum (this
### model is typically a bit less accurate, but much less complex, so
### better if you want to interpret the coefficients). mean(min): take
### the min of each K-CV error function, and then take their mean.
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
### example(ROChange).
 initial.regularization=0.001
### Passed to IntervalRegressionRegularized.
){
  validation.fold <- negative.auc <- threshold <- incorrect.labels <-
    variable <- value <- regularization <- folds <- status <- type <-
      vjust <- NULL
### The code above is to avoid CRAN NOTEs like
### IntervalRegressionCV: no visible binding for global variable
  n.observations <- check_features_targets(feature.mat, target.mat)
  stopifnot(is.integer(n.folds))
  stopifnot(is.integer(fold.vec))
  stopifnot(length(fold.vec) == n.observations)
  stopifnot(
    is.character(reg.type),
    length(reg.type)==1,
    reg.type %in% c("1sd", "mean(min)", "min(mean)"))
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
  validation.fold.vec <- unique(fold.vec)
  validation.data <- foreach(
    validation.fold=validation.fold.vec, .combine=rbind) %dopar% {
      ##print(validation.fold)
      is.validation <- fold.vec == validation.fold
      is.train <- !is.validation
      train.features <- feature.mat[is.train, all.finite, drop=FALSE]
      train.targets <- target.mat[is.train, , drop=FALSE]
      fit <- IntervalRegressionRegularized(
        train.features, train.targets, verbose=verbose,
        initial.regularization=initial.regularization)
      validation.features <- feature.mat[is.validation, , drop=FALSE]
      pred.log.lambda <- fit$predict(validation.features)
      validation.targets <- target.mat[is.validation, , drop=FALSE]
      too.small <- pred.log.lambda < validation.targets[, 1]
      too.big <- validation.targets[, 2] < pred.log.lambda
      is.error <- too.small | too.big
      left.term <- squared.hinge(pred.log.lambda-validation.targets[, 1])
      right.term <- squared.hinge(validation.targets[, 2]-pred.log.lambda)
      loss.vec <- colMeans(left.term+right.term)
      error.vec <- colSums(is.error)
      dt <- data.table(
        validation.fold,
        regularization=fit$regularization.vec,
        squared.hinge.loss=loss.vec,
        incorrect.intervals=error.vec)
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
    }
  variable.name <- if(!is.null(incorrect.labels.db)){
    "negative.auc"
  }else{
    "squared.hinge.loss"
  }
  vtall <- melt(
    validation.data,
    id.vars=c("validation.fold", "regularization"))
  variable.data <- vtall[variable==variable.name, ]
  stats <- variable.data[, list(
    mean=mean(value),
    sd=sd(value),
    folds=.N
    ), by=list(variable, regularization)][folds==max(folds),]
  min.each <- variable.data[, {
    .SD[which.min(value), ]
  }, by=validation.fold]
  min.mean <- stats[which.min(mean), ]
  upper.limit <- min.mean[, mean+sd]
  simplest.within.1sd <-
    stats[mean < upper.limit, ][which.max(regularization),]
  min.dt <- data.table(
    type=c("mean(min)", "min(mean)", "1sd"),
    vjust=c(1,2,1),
    regularization=c(
      mean(min.each$regularization),
      min.mean$regularization,
      simplest.within.1sd$regularization),
    variable=variable.name)
  min.dt[, status := ifelse(type == reg.type, "selected", "not")]
  fit <- IntervalRegressionRegularized(
    feature.mat, target.mat,
    initial.regularization=min.dt[status=="selected", regularization],
    factor.regularization=NULL,
    verbose=verbose)
  fit$plot.selectRegularization <- fit$plot <- ggplot()+
    ggtitle(paste0(
      "Regularization parameter selection using ",
      length(validation.fold.vec),
      "-fold cross-validation"
      ))+
    theme_bw()+
    geom_vline(aes(xintercept=-log(regularization)),
               data=min.dt[status=="selected",],
               color="grey",
               size=2)+
    geom_vline(aes(xintercept=-log(regularization), color=type),
               data=min.dt)+
    guides(color="none")+
    geom_text(aes(-log(regularization), max(variable.data$value),
                  vjust=vjust,
                  label=paste0(type, " "),
                  color=type),
              hjust=1,
              data=min.dt)+
    geom_segment(aes(
      -log(regularization), mean,
      xend=-log(min.mean$regularization), yend=mean,
      color=type),
               data=data.table(
                 simplest.within.1sd, type="1sd"))+
    theme_no_space()+
    facet_grid(variable ~ ., scales="free")+
    scale_color_manual(values=c(
                         "1sd"="red",
                         "mean(min)"="blue",
                         "min(mean)"="black"))+
    geom_ribbon(aes(
      -log(regularization),
      ymin=mean-sd,
      ymax=mean+sd),
                fill="grey",
                alpha=0.5,
                data=stats)+
    geom_line(aes(
      -log(regularization),
      mean,
      color="min(mean)"),
              data=stats)+
    geom_line(aes(-log(regularization), value, group=validation.fold),
              color="grey50",
              data=vtall[variable!="auc",])+
    geom_point(aes(
      -log(regularization),
      value,
      color="mean(min)"),
               data=min.each)+
    xlab("model complexity -log(regularization)")+
    ylab("")
  fit$plot.selectRegularization.data <- validation.data
  fit$plot.selectRegularization.vlines <- min.dt
  fit
}, ex=function(){

  if(interactive()){
    library(penaltyLearning)
    data("neuroblastomaProcessed", package="penaltyLearning", envir=environment())
    if(require(doParallel)){
      registerDoParallel()
    }
    set.seed(1)
    i.train <- 1:200
    fit <- with(neuroblastomaProcessed, IntervalRegressionCV(
      feature.mat[i.train,], target.mat[i.train,]))
    ## When only features and target matrices are specified for
    ## training, the squared hinge loss is used as the metric to
    ## minimize on the validation set.
    plot(fit)
    ## Create an incorrect labels data.table (first key is same as
    ## rownames of feature.mat and target.mat).
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
      "IntervalRegression model for regularization ",
      x$regularization.vec, 
      " with weights:\n",
      sep="")
    x <- t(x$pred.param.mat)
    rownames(x) <- ""
    print(x)
  }else{
    cat(
      "IntervalRegression models [",
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
### Error if the algorithm has not found an optimal solution after
### this many iterations.
 weight.vec=NULL,
### A numeric vector of weights for each training example.
 Lipschitz=NULL,
### A numeric scalar or NULL, which means to compute Lipschitz as the
### mean of the squared L2-norms of the rows of the feature matrix.
 verbose=2
### Cat messages: for restarts and at the end if >= 1, and for every
### iteration if >= 2.
 ){
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
  squared.hinge.deriv <- function(x){
    ifelse(x<1,2*(x-1),0)
  }  
  calc.loss <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge(linear.predictor-targets[,1])
    right.term <- squared.hinge(targets[,2]-linear.predictor)
    both.terms <- left.term+right.term
    weighted.loss.vec <- both.terms * weight.vec
    mean(weighted.loss.vec)
  }
  calc.grad <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge.deriv(linear.predictor-targets[,1])
    right.term <- squared.hinge.deriv(targets[,2]-linear.predictor)
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

