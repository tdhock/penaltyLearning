library(testthat)
context("ROChange")
library(penaltyLearning)
##library(data.table)

model.df <- structure(list(loss = c(11.9627022758106, 3.94540257494294, 2.94169197856002,
2.82924533907096, 2.72992938630457, 2.65399631185225, 2.55135727898647,
2.45833675556102, 2.37610825176777, 2.28285831363641, 2.2104078498086,
2.12689344888736, 2.05444298505956, 1.99401389635343, 1.92156343252562,
1.86727640796262, 1.81044827035527, 1.77872585064518, 1.72189771303783,
1.67261156766267), n.segments = 1:20), class = "data.frame", row.names = c(NA,
-20L))
modelSelection(model.df, complexity="n.segments")

data(neuroblastoma, package="neuroblastoma", envir=environment())
pid <- 81L
pro <- subset(neuroblastoma$profiles, profile.id==pid)
pro$pid <- pid
ann <- subset(neuroblastoma$annotations, profile.id==pid)
ann$pid <- pid
max.segments <- 20
segs.list <- list()
selection.list <- list()
chr <- "1"

for(chr in unique(ann$chromosome)){
  pro.chr <- subset(pro, chromosome==chr)
  print(chr)
  print(dim(pro.chr))
  fit <- jointseg::Fpsn(pro.chr$logratio, max.segments)
  str(fit)
  model.df <- data.frame(loss=fit$cost, n.segments=1:max.segments)
  print(model.df)
  dput(model.df)
  print(l <- modelSelectionC(model.df$loss, model.df$n.segments, 1:nrow(model.df)))
  selection.df <- modelSelection(model.df, complexity="n.segments")
  str(selection.df)
  selection.list[[chr]] <- data.frame(
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
    segs.list[[paste(chr, n.segments)]] <- data.frame(
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
  expect_equal(result$roc$fn, c(3, 2, 1, 0, 0, 0, 0))
  expect_equal(result$roc$fp, c(0, 0, 0, 0, 1, 2, 3))
})
test_that("perfect prediction achieved", {
  expect_equal(result$thresholds$errors, c(3, 0))
  expect_equal(result$auc, 1)
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

## Artificially tied penalty predictions. Although in real data sets
## it is possible to have the same predicted penalty (constant penalty
## function), it is extremely unlikely to have tied thresholds (that
## would mean the log.lambda values in modelSelection were the same
## for two data sets). Anyways, here is a test for that case.
some <- data.table(pid, chromosome=paste(c(1, 2)), n.segments=c(2, 1))
##error.list$model.errors[chromosome %in% c(1,2), .(chromosome, errors, n.segments, min.log.lambda)]
tie.pred <- data.table(pred)
tie.pred$pred.log.lambda[1:2] <- error.list$model.errors[some, on=list(
  pid, chromosome, n.segments), min.log.lambda - 2]
tie.result <- ROChange(error.list$model.errors, tie.pred, pvars)
test_that("six rows for six labels with one tie", {
  fn.vec <- c(3, 2, 1, 0, 0, 0)
  fp.vec <- c(0, 0, 0, 1, 2, 3)
  expect_equal(tie.result$roc$fp, fp.vec)
  expect_equal(tie.result$roc$fn, fn.vec)
  expect_equal(tie.result$roc$FPR, fp.vec/3)
  expect_equal(tie.result$roc$TPR, 1-fn.vec/3)
})
test_that("perfect prediction not achieved", {
  expect_equal(tie.result$thresholds$errors, c(3, 1))
  expect_equal(tie.result$auc, 17/18)
})

ann.breakBad <- data.table(ann)
ann.breakBad[, weight := ifelse(annotation=="breakpoint", 2, 1)]
error.breakBad <- labelError(
  selection, ann.breakBad, changes,
  problem.vars="chromosome", # for all three data sets.
  model.vars="n.segments", # for changes and selection.
  change.var="chromStart", # column of changes with breakpoint position.
  label.vars=c("min", "max")) # limit of labels in ann.
result.breakBad <- ROChange(error.breakBad$model.errors, tie.pred, "chromosome")
test_that("six rows for six labels with one tie, fn worse", {
  fn.vec <- c(6, 4, 2, 0, 0, 0)
  fp.vec <- c(0, 0, 0, 1, 2, 3)
  expect_equal(result.breakBad$roc$fp, fp.vec)
  expect_equal(result.breakBad$roc$fn, fn.vec)
  expect_equal(result.breakBad$roc$FPR, fp.vec/3)
  expect_equal(result.breakBad$roc$TPR, 1-fn.vec/6)
})
test_that("perfect prediction not achieved when fn worse", {
  expect_equal(result.breakBad$thresholds$errors, c(3, 1))
  expect_equal(result.breakBad$auc, 17/18)
  expect_equal(result.breakBad$thresholds$fp, c(3, 1))
  expect_equal(result.breakBad$thresholds$fn, c(0, 0))
})

ann.normBad <- data.table(ann)
ann.normBad[, weight := ifelse(annotation=="breakpoint", 1, 2)]
error.normBad <- labelError(
  selection, ann.normBad, changes,
  problem.vars="chromosome", # for all three data sets.
  model.vars="n.segments", # for changes and selection.
  change.var="chromStart", # column of changes with breakpoint position.
  label.vars=c("min", "max")) # limit of labels in ann.
result.normBad <- ROChange(error.normBad$model.errors, tie.pred, "chromosome")
test_that("six rows for six labels with one tie, fn worse", {
  fn.vec <- c(3, 2, 1, 0, 0, 0)
  fp.vec <- c(0, 0, 0, 2, 4, 6)
  expect_equal(result.normBad$roc$fp, fp.vec)
  expect_equal(result.normBad$roc$fn, fn.vec)
  expect_equal(result.normBad$roc$FPR, fp.vec/6)
  expect_equal(result.normBad$roc$TPR, 1-fn.vec/3)
})
test_that("perfect prediction not achieved when fn worse", {
  expect_equal(result.normBad$thresholds$errors, c(6, 1))
  expect_equal(result.normBad$auc, 17/18)
  expect_equal(result.normBad$thresholds$fp, c(6, 0))
  expect_equal(result.normBad$thresholds$fn, c(0, 1))
})

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

test_that("error for missing columns in model table", {
  expect_error({
    ROChange(data.table(chromosome="foo"), pred, "chromosome")
  }, "models should have columns")
})
