library(testthat)
context("IntervalRegression")
library(penaltyLearning)

test_that("IRCV converges on small data with prev problems", {
  data("notConverging", package="penaltyLearning")
  fit <- with(notConverging, IntervalRegressionCV(
    X.mat, y.mat, fold.vec=fold.vec,
    min.obs=nrow(X.mat),
    LAPPLY=lapply,
    verbose=0))
  expect_is(fit, "list")
})

feature.mat <- structure(c(
  1.79565978525711, 1.79606119914403, 1.71572069625467,
  1.61288792166534, 1.65564225297093, 1.66162385851159, 1.57115110896884,
  1.78610128168745), .Dim = c(8L, 1L), .Dimnames = list(c(
    "199_chr2",
    "479_chr2", "409_chr2", "18_chr3", "69_chr17", "68_chr4", "168_chr11",
    "551_chr3"), "log2.n"))
target.mat <- structure(c(
  -2.9286431049021, -2.98647204509594, -1.86813212772763,
  -Inf, -2.0471474145434, -4.31207055103816, -Inf, -2.12592617754227,
  Inf, Inf, Inf, 2.48942330405145, Inf, Inf, 1.58204507943406,
  Inf), .Dim = c(8L, 2L), .Dimnames = list(c(
    "199_chr2", "479_chr2",
    "409_chr2", "18_chr3", "69_chr17", "68_chr4", "168_chr11", "551_chr3"
  ), c("min.log.lambda", "max.log.lambda")))
fold.vec <- c(2L, 1L, 1L, 1L, 1L, 2L, 2L, 2L)
test_that("cv works for 8 observations", {
  fit <- IntervalRegressionCV(
    feature.mat, target.mat, fold.vec=fold.vec, min.obs=8)
  expect_is(fit, "list")
})

test_that("cv errors when 1sd specified but is not defined", {
  expect_error({
    IntervalRegressionCV(
      feature.mat, target.mat, fold.vec=fold.vec, min.obs=8, reg.type="1sd")
  }, "reg.type=1sd undefined; try another reg.type (min) or decrease initial.regularization",
  fixed=TRUE)
})

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

set.seed(1)
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
test_that("error when no upper limits", {
  expect_error({
    IntervalRegressionCV(
      train.feature.mat[no.upper,],
      train.target.mat[no.upper,],
      min.observations=sum(no.upper))
  }, "no upper limits")
})

fold.vec <- as.integer(ifelse(no.lower, 1, 2))
test_that("error when one fold has no upper/lower limits", {
  expect_error({
    IntervalRegressionCV(
      train.feature.mat,
      train.target.mat,
      fold.vec=fold.vec)
  }, "each fold should have at least one upper and one lower limit")
})
