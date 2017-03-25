library(testthat)
context("demo8")
library(penaltyLearning)

data(demo8, package="penaltyLearning")

set.seed(1)
fit <- with(demo8, IntervalRegressionCV(
  feature.mat, target.mat, min.observations=8))

test_that("CV model prints weights", {
  expect_output({
    print(fit)
  }, "IntervalRegression model for regularization")
})

test_that("CV model coef ncol 1", {
  expect_equal(ncol(coef(fit)), 1)
})

test_that("valid CV model for 8 train data with NA features", {
  pred.vec <- fit$predict(demo8$feature.mat)
  expect_equal(length(pred.vec), nrow(demo8$feature.mat))
  expect_true(is.numeric(pred.vec))
  expect_true(all(is.finite(pred.vec)))
})

test_that("plot(CV) returns ggplot", {
  gg <- plot(fit)
  expect_is(gg, "ggplot")
})  

fit <- with(demo8, IntervalRegressionRegularized(
  feature.mat, target.mat))
test_that("Regularized model contains plots and data", {
  expect_is(plot(fit), "ggplot")
  expect_is(fit$plot.residual, "ggplot")
  expect_is(fit$plot.weight, "ggplot")
  expect_is(fit$plot.residual.data, "data.table")
  expect_is(fit$plot.weight.data, "data.table")
})

test_that("Regularized model coef ncol", {
  expect_equal(ncol(coef(fit)), ncol(fit$pred.param.mat))
})

test_that("Regularized model prints summary", {
  expect_output({
    print(fit)
  }, "IntervalRegression models")
})

