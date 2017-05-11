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
set.seed(1)
fit <- IntervalRegressionCV(train.feature.mat, train.target.mat)
pred.vec <- predict(fit, neuroblastomaProcessed$feature.mat[allchr.rows, ])
roc.list <- targetIntervalROC(neuroblastomaProcessed$target.mat[allchr.rows, ], pred.vec)
test_that("perfect prediction for six profiles", {
  ##http://members.cbio.mines-paristech.fr/~thocking/change-tutorial/Supervised.html
  expect_equal(roc.list$auc, 1)
  expect_equal(roc.list$thresholds[threshold=="predicted", errors], 0)
})
