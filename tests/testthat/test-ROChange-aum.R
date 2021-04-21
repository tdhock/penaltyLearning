library(testthat)
context("ROChange")
library(penaltyLearning)
library(data.table)

ggcost <- function(){
  if(interactive() && require("ggplot2")){
    pred.dt <- data.table(predictions)
    gg.dt <- data.table(problem=pred.dt$problem)[, {
      data.table(log.pen=seq(-2, 2, by=0.5))[, {
        pred.dt$pred.log.lambda[problem] <- log.pen
        with(
          ROChange(models, pred.dt, "problem"),
          data.table(aum))
      }, by=log.pen]
    }, by=problem]
    ggplot()+
      geom_vline(aes(
        xintercept=pred.log.lambda),
        data=predictions)+
      geom_point(aes(
        log.pen, aum),
        data=gg.dt)+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_grid(problem ~ .)+
      coord_equal()
  }
}

models <- data.table(
  fp=c(1, 0, 1, 0, 0, 0,0,0),
  fn=c(0, 0, 0, 0, 0, 1, 0, 1),
  possible.fn=c(0,0,0,0,1,1,1,1),
  possible.fp=c(1,1,1,1,0,0,0,0),
  min.log.lambda=c(-Inf,-1, 0, 1,-Inf,-1,0,1),
  max.log.lambda=c(-1,0,1, Inf,-1,0,1,Inf),
  labels=1,
  problem=c(1,1,1,1,2,2,2,2))
models[, errors := fp+fn]
predictions <- data.table(problem=c(1,2), pred.log.lambda=c(1,0))
ggcost()
test_that("noncvx 1fp[-1,0] 1fn[0,1]", {
  L <- ROChange(models, predictions, "problem")
  ggplot()+geom_segment(aes(
    min.thresh, min.fp.fn, xend=max.thresh, yend=min.fp.fn),
    data=L$roc)
  expect_equal(L$aum, 1)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, c(1,2))
  expect_equal(L$aum.grad$lo, c(1,1))
  expect_equal(L$aum.grad$hi, c(-1,-1))
})

predictions <- data.table(problem=c(1,2), pred.log.lambda=c(0, Inf))
test_that("1fp[-1,0] 1fn[0,1]", {
  expect_error({
    ROChange(models, predictions, "problem")
  }, "all predictions must be finite")
})

models <- data.table(
  fp=c(1, 0, 0, 0),
  fn=c(0, 0, 0, 1),
  possible.fn=c(0,0,1,1),
  possible.fp=c(1,1,0,0),
  min.log.lambda=c(-Inf,0,-Inf,0),
  max.log.lambda=c(0,Inf,0,Inf),
  labels=1,
  problem=c(1,1,2,2))
models[, errors := fp+fn]
predictions <- data.table(problem=c(1,2), pred.log.lambda=0)
ggcost()
test_that("1fp[-1,0] 1fn[0,1]", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, c(1,2))
  expect_equal(L$aum.grad$lo, c(-1,0))
  expect_equal(L$aum.grad$hi, c(0,1))
})

predictions <- data.table(problem=c(1,2), pred.log.lambda=c(1, -1))
ggcost()
test_that("1fp[0,0] 1fn[0,0]", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, c(1,2))
  expect_equal(L$aum.grad$lo, c(0,0))
  expect_equal(L$aum.grad$hi, c(0,0))
})

models <- data.table(
  fp=c(1, 0, 0, 0, 0, 0),
  fn=c(0, 0, 0, 1, 0, 0),
  possible.fn=c(0,0,1,1,1,1),
  possible.fp=c(1,1,0,0,0,1),
  min.log.lambda=c(-Inf,0,-Inf,0,1, -Inf),
  max.log.lambda=c(0,Inf,0,1,Inf, Inf),
  labels=1,
  problem=c(1,1,2,2,2,3))
models[, errors := fp+fn]
predictions <- data.table(problem=c(1,2,3), pred.log.lambda=0)
ggcost()
test_that("1fp[-1,0] 1fn[0,1], no change", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, c(1,2,3))
  expect_equal(L$aum.grad$lo, c(-1,0,0))
  expect_equal(L$aum.grad$hi, c(0,1,0))
})

predictions <- data.table(problem=c(1,2), pred.log.lambda=c(1, -1))
ggcost()
test_that("three problems but two predictions", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, c(1,2))
  expect_equal(L$aum.grad$lo, c(0,0))
  expect_equal(L$aum.grad$hi, c(0,0))
})

predictions <- data.table(problem=c(1,2,2), pred.log.lambda=c(1, -1,0))
test_that("three models, three predictions with problem", {
  expect_error({
    ROChange(models, predictions, "problem")
  }, "more than one prediction per problem")
})

models <- data.table(
  fp=c(1, 0, 0, 0, 1, 0),
  fn=c(0, 0, 0, 2, 0, 0),
  possible.fn=c(0,0,2,2,0,0),
  possible.fp=c(1,1,0,0,1,1),
  min.log.lambda=c(-Inf,0,-Inf,0, -Inf,0),
  max.log.lambda=c(0,Inf,0,Inf,0,Inf),
  labels=c(1,1,2,2,1,1),
  problem=c(1,1,2,2,3,3))
models[, errors := fp+fn]
predictions <- data.table(problem=c(1,2,3), pred.log.lambda=0)
ggcost()
test_that("1fp[-1,0] 2fn[0,2] 1fp[-1,0]", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, 1:3)
  expect_equal(L$aum.grad$lo, c(-1,0,-1))
  expect_equal(L$aum.grad$hi, c(0,2,0))
})

models <- data.table(
  fp=c(1, 0, 0, 0, 1, 0),
  fn=c(0, 0, 0, 1, 0, 0),
  possible.fn=c(0,0,1,1,0,0),
  possible.fp=c(1,1,0,0,1,1),
  min.log.lambda=c(-Inf,0,-Inf,0, -Inf,0),
  max.log.lambda=c(0,Inf,0,Inf,0,Inf),
  labels=1,
  problem=c(1,1,2,2,3,3))
models[, errors := fp+fn]
predictions <- data.table(problem=c(1,2,3), pred.log.lambda=0)
ggcost()
test_that("1fp[-1,0] 1fn[0,1] 1fp[-1,0]", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, 1:3)
  expect_equal(L$aum.grad$lo, c(-1,0,-1))
  expect_equal(L$aum.grad$hi, c(0,1,0))
})

models <- data.table(
  fp=c(2, 0,   0, 0,   2, 1, 0, 0),
  fn=c(0, 0,   0, 1,   0, 0, 1, 2),
  possible.fn=c(0,0,1,1,2,2,2,2),
  possible.fp=c(2,2,0,0,2,2,2,2),
  min.log.lambda=c(-Inf,0,-Inf,0,-Inf, -1,0, 1),
  max.log.lambda=c(0,Inf,0,Inf,   -1, 0,1, Inf),
  labels=c(2,2,1,1,2,2,2,2),
  problem=c(1,1,2,2,3,3,3,3))
models[, errors := fp+fn]
predictions <- data.table(problem=c(1,2,3), pred.log.lambda=0)
ggcost()
test_that("2fp[-2,0] 1fn[0,1] 2fp2fn(0)[-1,1]", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 0)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, 1:3)
  expect_equal(L$aum.grad$lo, c(-2,0,-1))
  expect_equal(L$aum.grad$hi, c(0,1,1))
})

predictions <- data.table(problem=c(1,2,3), pred.log.lambda=c(0,0,1))
ggcost()
test_that("2fp[-2,-1] 1fn[0,1] 2fp2fn(1)[1,2]", {
  L <- ROChange(models, predictions, "problem")
  expect_equal(L$aum, 1)
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, 1:3)
  expect_equal(L$aum.grad$lo, c(-2,0,1))
  expect_equal(L$aum.grad$hi, c(-1,1,2))
})

models <- data.table(
  fp=c(4, 0,   0, 0,   2, 1, 0, 0),
  fn=c(0, 0,   0, 1,   0, 0, 1, 2),
  possible.fn=c(0,0,1,1,2,2,2,2),
  possible.fp=c(4,4,0,0,2,2,2,2),
  min.log.lambda=c(-Inf,0,-Inf,0,-Inf, -1,0, 1),
  max.log.lambda=c(0,Inf,0,Inf,   -1, 0,1, Inf),
  labels=c(4,4,1,1,2,2,2,2),
  problem=c(1,1,2,2,3,3,3,3))
models[, errors := fp+fn]
predictions <- data.table(
  problem=c(1,2,3), pred.log.lambda=c(-1,0,0))
ggcost()

test_that("4fp(-1)[-3,-2](0)[-2,0] 1fn[0,1] 2fp2fn(0)[-1,1](1)[1,2]", {
  L <- ROChange(models, predictions, "problem")
  print(L$aum.grad)
  expect_equal(L$aum.grad$problem, 1:3)
  expect_equal(L$aum.grad$lo, c(-3,1,1))
  expect_equal(L$aum.grad$hi, c(-2,1,2))
})

test_that("auc=2 for one error curve with one loop", {
  before.dt <- data.table(
    tp=0,
    fp=0,
    possible.tp=1,
    possible.fp=1)
  rep.dt <- data.table(
    tp=c(1, 1, 0, 0),
    fp=c(0, 1, 1, 0),
    possible.tp=1,
    possible.fp=1)
  after.dt <- data.table(
    tp=c(1, 1),
    fp=c(0, 1),
    possible.tp=1,
    possible.fp=1)
  rep.list <- replicate(1, rep.dt, simplify=FALSE)
  several.dt <- do.call(rbind, rep.list)
  segs.dt <- rbind(before.dt, several.dt, after.dt)
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
  expect_equal(L$auc, 2)
})

