library(testthat)
context("demo8")
library(penaltyLearning)

data(demo8, package="penaltyLearning")

test_that("valid CV model for 8 train data", {
  fit <- with(demo8, IntervalRegressionCV(
    feature.mat, target.mat, min.observations=8))
  pred.vec <- fit$predict(demo8$feature.mat)
  expect_equal(length(pred.vec), nrow(demo8$feature.mat))
})
