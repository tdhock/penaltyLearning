library(testthat)
context("features")
library(penaltyLearning)
data(neuroblastoma, package="neuroblastoma")

one <- subset(neuroblastoma$profiles, profile.id=="1" & chromosome=="1")
f.vec <- featureVector(one$logratio)

test_that("median absolute difference computed", {
  expect_equal(
    f.vec[["diff abs quantile.identity.50%"]],
    median(abs(diff(one$logratio))))
})

test_that("error for data.frame", {
  expect_error({
    featureVector(one)
  }, "data.vec must be a numeric data sequence with at least two elements, all of which are finite (not missing)", fixed=TRUE)
})

two <- subset(neuroblastoma$profiles, profile.id=="2" & chromosome=="2")
f2 <- featureVector(two$logratio)

test_that("feature vectors are the same size", {
  expect_equal(length(f2), length(f.vec))
})

three <- subset(neuroblastoma$profiles, profile.id %in% 1:3)
f.mat <- featureMatrix(three, c("profile.id", "chromosome"), "logratio")

test_that("feature matrix has same columns as vector", {
  expect_identical(colnames(f.mat), names(f2))
})

u3 <- with(three, unique(paste(profile.id, chromosome)))
test_that("feature matrix has expected row names", {
  expect_identical(rownames(f.mat), u3)
})
