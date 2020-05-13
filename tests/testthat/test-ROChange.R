library(testthat)
context("ROChange")
library(penaltyLearning)
library(data.table)

test_that("AUC of reverse ROC curve is 1", {
  segs.dt <- data.table(
    tp=c(0, 1, 1),
    fp=c(0, 0, 1),
    possible.tp=1,
    possible.fp=1)
  n.breaks <- nrow(segs.dt)-1L
  break.vec <- 1:n.breaks
  segs.dt[, min.log.lambda := c(-Inf, break.vec)]
  segs.dt[, max.log.lambda := c(break.vec, Inf)]
  segs.dt[, problem := 1]
  segs.dt[, fn := possible.tp-tp]
  segs.dt[, possible.fn := possible.tp]
  segs.dt[, errors := fp+fn]
  segs.dt[, labels := 1]
  pred.dt <- data.table(pred.log.lambda=1.5, problem=1)
  L <- ROChange(segs.dt, pred.dt, "problem")
  expect_equal(L$auc, 1)
})

test_that("error for labels less than errors", {
  segs.dt <- data.table(
    tp=c(1, 2, 2),
    fp=c(1, 1, 2),
    possible.tp=2,
    possible.fp=2)
  n.breaks <- nrow(segs.dt)-1L
  break.vec <- 1:n.breaks
  segs.dt[, min.log.lambda := c(-Inf, break.vec)]
  segs.dt[, max.log.lambda := c(break.vec, Inf)]
  segs.dt[, problem := 1]
  segs.dt[, fn := possible.tp-tp]
  segs.dt[, possible.fn := possible.tp]
  segs.dt[, errors := fp+fn]
  segs.dt[, labels := 1]
  pred.dt <- data.table(pred.log.lambda=1.5, problem=1)
  expect_error({
    ROChange(segs.dt, pred.dt, "problem")
  }, "errors should be in [0,labels]", fixed=TRUE)
})

test_that("AUC of reverse incomplete ROC curve is 1", {
  segs.dt <- data.table(
    tp=c(1, 2, 2),
    fp=c(0, 0, 1),
    possible.tp=2,
    possible.fp=2)
  n.breaks <- nrow(segs.dt)-1L
  break.vec <- 1:n.breaks
  segs.dt[, min.log.lambda := c(-Inf, break.vec)]
  segs.dt[, max.log.lambda := c(break.vec, Inf)]
  segs.dt[, problem := 1]
  segs.dt[, fn := possible.tp-tp]
  segs.dt[, possible.fn := possible.tp]
  segs.dt[, errors := fp+fn]
  segs.dt[, labels := 2]
  pred.dt <- data.table(pred.log.lambda=1.5, problem=1)
  L <- ROChange(segs.dt, pred.dt, "problem")
  L$roc
  expect_equal(L$auc, 1)
})

simple.err <- rbind(
  data.table(
    pid=81, chromosome="1", min.log.lambda=c(-Inf, 0), max.log.lambda=c(0, Inf),
    errors=c(0, 1), labels=1,
    fn=c(0, 1), possible.fn=1,
    fp=0, possible.fp=0),
  data.table(
    pid=81, chromosome="2", min.log.lambda=c(-Inf, 1), max.log.lambda=c(1, Inf),
    errors=c(1, 0), labels=1,
    fn=0, possible.fn=0,
    fp=c(1, 0), possible.fp=1))
test_that("only one prediction row even when prediction is on threshold", {
  L <- ROChange(simple.err, ok.pred, pvars)
  pred.dt <- L$thresholds[threshold=="predicted"]
  expect_equal(nrow(pred.dt), 1)
})

two.pred <- rbind(ok.pred, ok.pred)
test_that("two predictions for the same problem is an error", {
  expect_error({
    ROChange(error.list$model.errors, two.pred, pvars)
  }, "more than one prediction per problem")
})

test_that("error for missing columns in model table", {
  expect_error({
    ROChange(data.table(chromosome="foo"), pred, "chromosome")
  }, "models should have columns")
})

