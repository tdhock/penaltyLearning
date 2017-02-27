library(testthat)
context("errors")
library(penaltyLearning)

data(neuroblastomaProcessed, package="penaltyLearning")
fit <- with(neuroblastomaProcessed, IntervalRegressionUnregularized(
  feature.mat[, c("log.hall", "log.n")], target.mat))
no.name.mat <- matrix(c(1, 1, 1, NA), 2, 2)
test_that("predict without all pred col names is an error", {
  expect_error({
    fit$predict(no.name.mat)
  }, "need some missing features for prediction: log.hall, log.n")
})

named.mat <- no.name.mat
colnames(named.mat) <- c("log.hall", "log.n")
test_that("predict NA results in NA", {
  pred.vec <- fit$predict(named.mat)
  expect_identical(as.logical(is.na(pred.vec)), c(FALSE, TRUE))
})
