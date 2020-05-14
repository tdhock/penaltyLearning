library(testthat)
context("ROChange-no-thresh")
library(penaltyLearning)
library(data.table)

m <- function(problem, min.log.lambda, max.log.lambda, errors, fp, fn, labels, possible.fp, possible.fn){
  data.table(problem, min.log.lambda, max.log.lambda, errors, fp, fn, labels, possible.fp, possible.fn)
}
model.dt <- rbind(#           Er fp fn N  fp fn
  m("two-thresh", -Inf, -500, 1, 1, 0, 1, 1, 1),
  m("two-thresh", -500, 500,  0, 0, 0, 1, 1, 1),
  m("two-thresh", 500, Inf,   1, 0, 1, 1, 1, 1),
  m("always-fp", -Inf, Inf,   0, 1, 0, 1, 1, 0),
  m("always-fn", -Inf, Inf,   0, 0, 1, 1, 0, 1),
  m("always-tn", -Inf, Inf,   0, 0, 0, 1, 1, 0),
  m("always-tp", -Inf, Inf,   0, 0, 0, 1, 0, 1))
pred.dt <- data.table(
  problem=unique(model.dt$problem),
  pred.log.lambda=0)
test_that("problem with no thresh is OK", {
  L <- ROChange(model.dt, pred.dt, "problem")
  expect_is(L, "list")
  expect_equal(L$roc$fp, fp <- c(2, 1, 1))
  expect_equal(L$roc$fn, fn <- c(1, 1, 2))
  expect_equal(L$roc$FPR, fp/3)
  expect_equal(L$roc$TPR, 1-fn/3)
})

test_that("inconsistent possible.fn/possible.fp/labels is an error", {
  for(col.name in c("possible.fn", "possible.fp", "labels")){
    inconsistent.dt <- data.table(model.dt)
    inconsistent.dt[[col.name]][1] <- inconsistent.dt[[col.name]][1]+1
    msg <- paste(
      col.name,
      "should be constant for each problem")
    expect_error({
      ROChange(inconsistent.dt, pred.dt, "problem")
    }, msg)
  }
})

test_that("negative possible.fn/possible.fp/labels is an error", {
  for(col.name in c("possible.fn", "possible.fp", "labels")){
    negative.dt <- data.table(model.dt)
    negative.dt[[col.name]][1:3] <- -1
    expect_error({
      ROChange(negative.dt, pred.dt, "problem")
    }, "possible.fn/possible.fp/labels should be non-negative")
  }
})

test_that("missing data is an error", {
  for(col.name in names(model.dt)){
    missing.dt <- data.table(model.dt)
    missing.dt[[col.name]][1] <- NA
    msg <- paste(col.name, "should not be NA")
    expect_error({
      ROChange(missing.dt, pred.dt, "problem")
    }, msg)
  }
})

test_that("problem with no Inf max.log.lambda is an error", {
  no.Inf <- data.table(model.dt)
  no.Inf[.N, max.log.lambda := -600]
  expect_error({
    ROChange(no.Inf, pred.dt, "problem")
  }, "for every problem, the smallest min.log.lambda should be -Inf, and the largest max.log.lambda should be Inf")
})

test_that("problem with no -Inf min.log.lambda is an error", {
  no.Inf <- data.table(model.dt)
  no.Inf[1, min.log.lambda := -600]
  expect_error({
    ROChange(no.Inf, pred.dt, "problem")
  }, "for every problem, the smallest min.log.lambda should be -Inf, and the largest max.log.lambda should be Inf")
})

test_that("problem with inconsistent min/max.log.lambda is an error", {
  inconsistent.dt <- data.table(model.dt)
  inconsistent.dt[1, max.log.lambda := -6]
  expect_error({
    ROChange(inconsistent.dt, pred.dt, "problem")
  }, "max.log.lambda should be equal to the next min.log.lambda")
})

possible.name.vec <- c(
  errors="labels",
  fp="possible.fp",
  fn="possible.fn")
test_that("problem with fp/fn/errors out of range is an error", {
  for(col.name in c("fn", "fp", "errors")){
    out.dt <- data.table(model.dt)
    for(bad.value in c(-1, 2)){
      out.dt[[col.name]][1] <- bad.value
      poss.name <- possible.name.vec[[col.name]]
      msg <- paste0(
        col.name,
        " should be in [0,",
        poss.name,
        "]")
      expect_error({
        ROChange(out.dt, pred.dt, "problem")
      }, msg, fixed=TRUE)
    }
  }
})



