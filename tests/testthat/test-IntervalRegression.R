library(testthat)
context("IntervalRegression")
library(penaltyLearning)

followup <- data.frame(
  profile.id=c(10L, 8L, 4L, 6L, 11L, 1L),
  status=c("ok", "relapse", "relapse", "ok", "ok", "relapse"))
data(neuroblastomaProcessed)
allchr.pattern <- paste0("^(", paste(followup$profile.id, collapse="|"), ")[.]")
allchr.rows <- grep(allchr.pattern, rownames(neuroblastomaProcessed$target.mat), value=TRUE)
train.rows <- grep("[.]11$", allchr.rows, invert=TRUE, value=TRUE)
train.feature.mat <- neuroblastomaProcessed$feature.mat[train.rows, ]
train.target.mat <- neuroblastomaProcessed$target.mat[train.rows, ]
set.seed(2)
fit <- IntervalRegressionCV(train.feature.mat, train.target.mat)
pred.vec <- predict(fit, neuroblastomaProcessed$feature.mat[allchr.rows, ])
roc.list <- targetIntervalROC(neuroblastomaProcessed$target.mat[allchr.rows, ], pred.vec)
test_that("perfect prediction for six profiles", {
  ##http://members.cbio.mines-paristech.fr/~thocking/change-tutorial/Supervised.html
  expect_equal(roc.list$auc, 1)
  expect_equal(roc.list$thresholds[threshold=="predicted", errors], 0)
})

test_that("threshold argument passed to IntervalRegressionInternal", {
  fit <- IntervalRegressionCV(
    train.feature.mat, train.target.mat, threshold=0.1)
})

no.lower <- train.target.mat[,1] == -Inf
test_that("error when no lower limits", {
  expect_error({
    IntervalRegressionCV(
      train.feature.mat[no.lower,],
      train.target.mat[no.lower,],
      min.observations=sum(no.lower))
  }, "no lower limits")
})

no.upper <- train.target.mat[,2] == Inf
test_that("error when no lower limits", {
  expect_error({
    IntervalRegressionCV(
      train.feature.mat[no.upper,],
      train.target.mat[no.upper,],
      min.observations=sum(no.upper))
  }, "no upper limits")
})
