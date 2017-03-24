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
  }, "columns needed for prediction but not present: log.hall, log.n")
})

named.mat <- no.name.mat
colnames(named.mat) <- c("log.hall", "log.n")
test_that("predict NA results in NA", {
  pred.vec <- predict(fit, named.mat)#also tests S3 method.
  expect_identical(as.logical(is.na(pred.vec)), c(FALSE, TRUE))
})

test_that("S3 print method returns numeric matrix with dimnames", {
  expect_is(fit, "IntervalRegression")
  result <- print(fit)
  expect_is(result, "matrix")
  expect_true(is.numeric(result))
  expect_identical(dimnames(result), list(
    c("(Intercept)", "log.hall", "log.n"),
    "0"))
})
