library(testthat)
context("ROChange-jointseg")
library(penaltyLearning)
library(data.table)

if(requireNamespace("neuroblastoma") && requireNamespace("jointseg")){
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  ## Find some profiles with all negative or positive labels.
  all.ann <- data.table(neuroblastoma$annotations)
  cast.ann <- dcast(
    all.ann, profile.id ~ annotation, value.var="annotation", fun.aggregate=length)
  ## profile.id 4 has all breakpoint labels.
  pid <- 4
  pro <- subset(neuroblastoma$profiles, profile.id==pid)
  ann <- all.ann[profile.id==pid, ]
  max.segments <- 20
  segs.list <- list()
  selection.list <- list()
  for(chr in unique(ann$chromosome)){
    pro.chr <- subset(pro, chromosome==chr)
    fit <- jointseg::Fpsn(
      pro.chr$logratio, max.segments)
    model.df <- data.frame(loss=fit$J.est, n.segments=1:max.segments)
    selection.df <- modelSelection(model.df, complexity="n.segments")
    selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
    for(n.segments in 1:max.segments){
      end <- fit$t.est[n.segments, 1:n.segments]
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
        chromEnd)
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
  ## The BIC model selection criterion is lambda = log(n), where n is
  ## the number of data points to segment. This implies log(lambda) =
  ## log(log(n)) = the log2.n feature in all.features.mat.
  pred <- pro.with.ann[, list(pred.log.lambda=log(log(.N))), by=chromosome]
  test_that("error when no negative labels", {
    expect_error({
      ROChange(error.list$model.errors, pred, "chromosome")
    }, "no negative labels")
  })
  ## profile.id 6 has no normal labels.
  pid <- 6
  pro <- subset(neuroblastoma$profiles, profile.id==pid)
  ann <- all.ann[profile.id==pid, ]
  max.segments <- 20
  segs.list <- list()
  selection.list <- list()
  for(chr in unique(ann$chromosome)){
    pro.chr <- subset(pro, chromosome==chr)
    fit <- jointseg::Fpsn(
      pro.chr$logratio, max.segments)
    model.df <- data.frame(loss=fit$J.est, n.segments=1:max.segments)
    selection.df <- modelSelection(model.df, complexity="n.segments")
    selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
    for(n.segments in 1:max.segments){
      end <- fit$t.est[n.segments, 1:n.segments]
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
        chromEnd)
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
  ## The BIC model selection criterion is lambda = log(n), where n is
  ## the number of data points to segment. This implies log(lambda) =
  ## log(log(n)) = the log2.n feature in all.features.mat.
  pred <- pro.with.ann[, list(pred.log.lambda=log(log(.N))), by=chromosome]
  test_that("error when no positive labels", {
    expect_error({
      ROChange(error.list$model.errors, pred, "chromosome")
    }, "no positive labels")
  })
  pid <- 81L
  pro <- subset(neuroblastoma$profiles, profile.id==pid)
  pro$pid <- pid
  ann <- subset(neuroblastoma$annotations, profile.id==pid)
  ann$pid <- pid
  segs.list <- list()
  selection.list <- list()
  for(chr in unique(pro$chr)){
    pro.chr <- subset(pro, chromosome==chr)
    max.segments <- min(20, nrow(pro.chr))
    fit <- jointseg::Fpsn(
      pro.chr$logratio, max.segments)
    model.df <- data.frame(loss=fit$J.est, n.segments=1:max.segments)
    selection.df <- modelSelection(model.df, complexity="n.segments")
    selection.list[[chr]] <- data.table(
      pid,
      chromosome=chr, selection.df)
    for(n.segments in 1:max.segments){
      cat(sprintf("chr=%s segments=%d\n", chr, n.segments))
      end <- fit$t.est[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
      (pro.chr$position[data.before.change]+
         pro.chr$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro.chr$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro.chr$position))
      segs.list[[paste(chr, n.segments)]] <- data.table(
        pid,
        chromosome=chr,
        n.segments,
        start,
        end,
        chromStart,
        chromEnd)
    }
  }
  segs <- do.call(rbind, segs.list)
  selection <- do.call(rbind, selection.list)
  changes <- segs[1 < start,]
  pvars <- c("chromosome", "pid")
  error.list <- labelError(
    selection, ann, changes,
    problem.vars=pvars, # for all three data sets.
    model.vars="n.segments", # for changes and selection.
    change.var="chromStart", # column of changes with breakpoint position.
    label.vars=c("min", "max")) # limit of labels in ann.
  pro.with.ann <- data.table(pro)[chromosome %in% ann$chromosome, ]
  bad.pred <- pro.with.ann[, list(
    pred.log.penalty=log(log(.N))
  ), by=pvars]
  test_that("informative error for no pred.log.lambda column", {
    expect_error({
      ROChange(error.list$model.errors, bad.pred, pvars)
    }, "predictions should be a data.frame with at least one row and a column named pred.log.lambda")
  })
  ## The BIC model selection criterion is lambda = log(n), where n is
  ## the number of data points to segment. This implies log(lambda) =
  ## log(log(n)) = the log2.n feature in all.features.mat.
  pred <- pro.with.ann[, list(
    pred.log.lambda=log(log(.N))
  ), by=pvars]
  result <- ROChange(error.list$model.errors, pred, pvars)
  test_that("seven rows for six labels", {
    expect_equal(result$roc$fp, c(3, 2, 1, 0, 0, 0, 0))
    expect_equal(result$roc$fn, c(0, 0, 0, 0, 1, 2, 3))
  })
  test_that("(FPR=1, TPR=0) in polygon but not roc", {
    expect_equal(result$auc.polygon[, sum(FPR==1 & TPR==0)], 1)
    expect_equal(result$roc[, sum(FPR==1 & TPR==0)], 0)
  })
  bad.pred <- data.table(
    pid,
    chromosome=paste(1:24),
    pred.log.lambda=0)
  error.list$model.errors[, table(chromosome)]
  test_that("informative error when predicting for unlabeled data", {
    expect_error({
      ROChange(error.list$model.errors, bad.pred, pvars)
    }, "some predictions do not exist in models")
  })

  ok.pred <- bad.pred[chromosome %in% c("1", "2")]
  test_that("predicting for one positive and one negative label is OK", {
    L <- ROChange(error.list$model.errors, ok.pred, pvars)
    expect_is(L, "list")
    expect_equal(L$thresholds$possible.fp, c(1, 1))
    expect_equal(L$thresholds$possible.fn, c(1, 1))
  })
}

