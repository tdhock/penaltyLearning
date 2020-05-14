featureVector <- structure(function
### Compute a feature vector of constant length which can be used as
### an input for supervised penalty learning. The output is a target
### interval of log(penalty) values that achieve minimum incorrect
### labels (see targetIntervals).
(data.vec
### numeric vector of ordered data.
){
  if(!(
    is.numeric(data.vec) &&
    1 < length(data.vec) &&
    all(is.finite(data.vec))
  )){
    stop("data.vec must be a numeric data sequence with at least two elements, all of which are finite (not missing)")
  }
  mu <- mean(data.vec)
  residual.vec <- data.vec - mu
  n <- length(data.vec)
  vec.list <- list(
    data=data.vec,
    residual=residual.vec,
    diff=diff(data.vec))
  nonlinear <- function(x){
    suppressWarnings(c(
      identity=x,
      sqrt=sqrt(x),
      log=log(x),
      loglog=log(log(x)),
      square=x*x))
  }
  feature.list <- list(n=nonlinear(n))
  for(data.type in names(vec.list)){
    vec <- vec.list[[data.type]]
    trans.list <- list(
      orig=vec,
      abs=abs(vec),
      square=vec * vec)
    for(trans.name in names(trans.list)){
      trans.vec <- trans.list[[trans.name]]
      feature.list[[paste(data.type, trans.name)]] <- nonlinear(c(
        sum=sum(trans.vec),
        mean=mean(trans.vec),
        sd=sd(trans.vec),
        quantile=quantile(trans.vec)))
    }
  }
  do.call(c, feature.list)
### Numeric vector of features.
}, ex=function(){

  x <- rnorm(10)
  penaltyLearning::featureVector(x)
  if(requireNamespace("neuroblastoma")){
    data(neuroblastoma, package="neuroblastoma", envir=environment())
    one <- subset(neuroblastoma$profiles, profile.id=="1" & chromosome=="1")
    (f.vec <- penaltyLearning::featureVector(one$logratio))
  }

})

featureMatrix <- structure(function
### Compute a feature matrix (segmentation problems x features).
(data.sequences,
### data.frame of sorted sequences of data to segment.
  problem.vars,
### character vector of columns of data.sequences to treat as
### segmentation problem IDs.
  data.var
### character vector of length 1 (column of data.sequences to treat as
### data to segment).
){
  if(!(
    is.data.frame(data.sequences) &&
    1 < nrow(data.sequences)
  )){
    stop("data.sequences must be a data.frame with at least two rows")
  }
  if(!(
    is.character(problem.vars) &&
      0 < length(problem.vars) &&
      all(!is.na(problem.vars)) &&
      all(problem.vars %in% names(data.sequences))
  )){
    stop("problem.vars must be a character vector of column names of data.sequences (IDs for separate segmentation problems)")
  }
  if(!(
    is.character(data.var) &&
    1 == length(data.var) &&
    !is.na(data.var) &&
    data.var %in% names(data.sequences)
  )){
    stop("data.var must be a character vector of length 1 (column name of data.sequences to treat as data to segment)")
  }
  feature.dt <- data.table(data.sequences)[, list(
    features=list(featureVector(.SD[[data.var]]))
  ), by=problem.vars]
  f.mat <- do.call(rbind, feature.dt$features)
  prob.mat <- as.matrix(feature.dt[, problem.vars, with=FALSE])
  rownames(f.mat) <- apply(prob.mat, 1, paste, collapse=" ")
  f.mat
### Numeric feature matrix. Some entries may be missing or infinite;
### these columns should be removed before model training.
}, ex=function(){

  test.df <- data.frame(
    id=rep(1:2, each=10),
    x=rnorm(20))
  penaltyLearning::featureMatrix(test.df, "id", "x")
  if(requireNamespace("neuroblastoma")){
    data(neuroblastoma, package="neuroblastoma", envir=environment())
    one <- subset(neuroblastoma$profiles, profile.id %in% c(1,2))
    f.mat <- penaltyLearning::featureMatrix(
      one, c("profile.id", "chromosome"), "logratio")
  }

})
