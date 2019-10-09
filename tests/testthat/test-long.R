library(testthat)
library(penaltyLearning)
context('long')
data("neuroblastomaProcessed", package="penaltyLearning")
X.mat <- neuroblastomaProcessed$feature.mat
y.mat <- neuroblastomaProcessed$target.mat
test_that("error for un-logged outputs", {
  expect_error({
    IntervalRegressionCV(X.mat, exp(y.mat), verbose=1, LAPPLY=lapply)
  }, "all targets are non-negative, and there is a big change between quantiles, so outputs are probably un-logged; this function expects real-valued (possibly negative) limits, so please try taking the log of your target matrix, or using check.unlogged=FALSE", fixed=TRUE)
})
