library(testthat)
context("modelSelection")
library(penaltyLearning)
data(oneSkip)

test_that("output intervals computed correctly", {
  df <- with(oneSkip$input, modelSelectionC(error, segments, peaks))
  expect_identical(df$model.complexity, oneSkip$output$model.complexity)
  expect_identical(df$min.lambda, oneSkip$output$min.lambda)
})

library(neuroblastoma)
library(Segmentor3IsBack)
data(neuroblastoma)
one <- subset(neuroblastoma$profiles, profile.id==599 & chromosome=="14")
max.segments <- 1000
fit <- Segmentor(one$logratio, model=2, Kmax=max.segments)
lik.df <- data.frame(lik=fit@likelihood, segments=1:max.segments)
pathR <- with(lik.df, modelSelectionR(lik, segments, segments))
pathC <- with(lik.df, modelSelectionC(lik, segments, segments))
test_that("C code agrees with R code for big data set", {
  expect_identical(pathR, pathC)
})

