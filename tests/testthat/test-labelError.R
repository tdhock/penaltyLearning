library(testthat)
context("labelError")
library(penaltyLearning)
library(data.table)

## Trivial edge cases.
## 123456789
##  (-0-]
##      (1]
##  | OK
##   | FP
##      | FP
##       | TP
##        | TP
##         | FN
label <- function(annotation, start, end){
  data.frame(prob="five", start, end, annotation)
}
ann.trivial <- rbind(
  label("1change", 6, 8),
  label("0changes", 2, 6))
models <- data.table(
  prob="five",
  algo="opart",
  complexity=c(-1, -3, -5, -6))
changes <- data.table(
  prob="five",
  pos=c(1, 7, 1, 6, 17, 11),
  algo="opart",
  complexity=c(-3, -5, -5, -6, -6, -6))
test_that("labelError throws informative errors", {
  expect_error({
    labelError(models, ann.trivial, changes)
  }, "problem.vars should be a character vector of column names present in models, changes, and labels (ID for separate changepoint detection problems)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob")
  }, "label.vars should be a 2-element character vector of labels column names (start and end of labeled region)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=character())
  }, "label.vars should be a 2-element character vector of labels column names (start and end of labeled region)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("foo", "bar"))
  }, "label.vars should be a 2-element character vector of labels column names (start and end of labeled region)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"))
  }, "change.var should be a column name of changes (position of predicted changepoints)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var=c("foo1"))
  }, "change.var should be a column name of changes (position of predicted changepoints)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var=c("pos", "end"))
  }, "change.var should be a column name of changes (position of predicted changepoints)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var="pos")
  }, "model.vars should be a column name of both models and changes (ID for model complexity, typically the number of changepoints or segments)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var="pos")
  }, "model.vars should be a column name of both models and changes (ID for model complexity, typically the number of changepoints or segments)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("end", "start"),
               change.var="pos",
               model.vars="complexity")
  }, "label start must be less than end", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start","end"),
               change.var="pos",
               model.vars=c("foo","complexity"))
  }, "model.vars should be a column name of both models and changes (ID for model complexity, typically the number of changepoints or segments)", fixed=TRUE)
})

trivial.list <- labelError(
  models, ann.trivial, changes,
  problem.vars="prob",
  label.vars=c("start", "end"),
  change.var="pos",
  model.vars=c("algo","complexity"))
test_that("1 TP for complexity=-5, 2 errors for -6", {
  trivial.list$model.errors[, {
    expect_equal(complexity, c(-1, -3, -5, -6))
    expect_equal(labels, c(2, 2, 2, 2))
    expect_equal(errors, c(1, 1, 0, 2))
    expect_equal(possible.fp, c(2, 2, 2, 2))
    expect_equal(fp, c(0, 0, 0, 1))
    expect_equal(possible.fn, c(1, 1, 1, 1))
    expect_equal(fn, c(1, 1, 0, 1))
  }]
})

ann.overlap <- rbind(
  label("1change", 5, 8),
  label("0changes", 2, 6))
test_that("error for overlapping labels", {
  expect_error({
    labelError(
      models, ann.overlap, changes,
      problem.vars="prob",
      label.vars=c("start", "end"),
      change.var="pos",
      model.vars="complexity")
  }, "each label end must be <= next label start", fixed=TRUE)
})

ann.unrecognized <- rbind(
  label("oneChange", 5, 8),
  label("0changes", 2, 4))
test_that("error for unrecognized labels", {
  expect_error({
    labelError(
      models, ann.unrecognized, changes,
      problem.vars="prob",
      label.vars=c("start", "end"),
      change.var="pos",
      model.vars="complexity")
  }, "labels$annotation must be one of annotations$annotation", fixed=TRUE)
})

test_that("label error works when some model cols are NA", {
  model.dt <- data.table(prob="five", model="foo", bar=NA)
  change.dt <- data.table(prob=character(), model=character(), pos=integer())
  label.dt <- label("1change", 5, 8)
  err.list <- labelError(
    model.dt, label.dt, change.dt,
    problem.vars="prob",
    label.vars=c('start', 'end'),
    change.var="pos",
    model.vars="model")
  expect_equal(err.list$model.errors$errors, 1)
})
