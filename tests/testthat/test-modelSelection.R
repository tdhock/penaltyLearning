library(testthat)
context("modelSelection")
library(penaltyLearning)
data(oneSkip)

test_that("output intervals computed correctly", {
  df <- with(oneSkip$input, modelSelection(error, segments, peaks))
  expect_identical(df$model.complexity, oneSkip$output$model.complexity)
})
