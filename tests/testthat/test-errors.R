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
pred.vec <- predict(fit, named.mat)#also tests S3 method.
test_that("predict NA results in NA", {
  expect_identical(as.logical(is.na(pred.vec)), c(FALSE, TRUE))
})

ok.targets <- rbind(
  c(1, 2),
  c(3, 4))
test_that("error for missing prediction", {
  expect_error({
    targetIntervalROC(ok.targets, pred.vec)
  }, "missing")
})

test_that("error for infinite prediction", {
  expect_error({
    targetIntervalROC(ok.targets, c(1, Inf))
  }, "infinite")
})

test_that("error for wrong size prediction", {
  expect_error({
    targetIntervalROC(ok.targets, 1:3)
  }, "length(pred) must be same as nrow(target.mat)", fixed=TRUE)
})

test_that("error for no finite limit", {
  expect_error({
    targetIntervalROC(rbind(ok.targets, c(-Inf, Inf)), 1:3)
  }, "each row of target.mat must have at least one finite limit")
})

ok3 <- rbind(ok.targets, c(-1, 5))
test_that("no error for finite limits", {
  roc.list <- targetIntervalROC(ok3, 1:3)
  expect_is(roc.list$thresholds, "data.table")
  expect_is(roc.list$auc.polygon, "data.table")
  expect_is(roc.list$roc, "data.table")
  expect_is(roc.list$auc, "numeric")
})

test_that("S3 coef method returns numeric matrix with dimnames", {
  expect_is(fit, "IntervalRegression")
  result <- coef(fit)
  expect_is(result, "matrix")
  expect_true(is.numeric(result))
  expect_identical(dimnames(result), list(
    c("(Intercept)", "log.hall", "log.n"),
    "0"))
})

test_that("error for more than two columns in target Unregularized", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionUnregularized(
      target.mat, feature.mat))
  }, "target.mat should be a numeric matrix with two columns (lower and upper limits of correct outputs)", fixed=TRUE)
})

test_that("error for more than two columns in target CV", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionCV(
      target.mat, feature.mat))
  }, "target.mat should be a numeric matrix with two columns (lower and upper limits of correct outputs)", fixed=TRUE)
})

test_that("error for different number of rows", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionCV(
      feature.mat[1:100,], target.mat[1:200,]))
  }, "feature.mat and target.mat should have the same number of rows", fixed=TRUE)
})

test_that("error for constant features", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionCV(
      cbind(x=rep(1, 200)), target.mat[1:200,]))
  }, "after filtering NA and constant features, none remain for training")
})

test_that("error for NA features", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionCV(
      cbind(x=rep(NA_real_, 200)), target.mat[1:200,]))
  }, "after filtering NA features, none remain for training")
})

test_that("error for NA and constant features", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionCV(
      cbind("NA"=rep(NA_real_, 200), constant=rep(1, 200)), target.mat[1:200,]))
  }, "after filtering NA and constant features, none remain for training")
})

test_that("error for un-named features", {
  expect_error({
    with(neuroblastomaProcessed, IntervalRegressionCV(
      cbind(rnorm(200)), target.mat[1:200,]))
  }, "feature.mat should be a numeric matrix with colnames (input features)", fixed=TRUE)
})
