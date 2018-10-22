library(testthat)
context("ROChange-no-thresh")
library(penaltyLearning)

m <- function(problem, min.log.lambda, max.log.lambda, errors, fp, fn, labels, possible.fp, possible.fn){
  data.table(problem, min.log.lambda, max.log.lambda, errors, fp, fn, labels, possible.fp, possible.fn)
}
model.dt <- rbind(#           Er fp fn N  fp fn 
  m("always-fp", -Inf, Inf,   0, 1, 0, 1, 1, 0),
  m("always-fn", -Inf, Inf,   0, 0, 1, 1, 0, 1),
  m("always-tn", -Inf, Inf,   0, 0, 0, 1, 1, 0),
  m("always-tp", -Inf, Inf,   0, 0, 0, 1, 0, 1),
  m("two-thresh", -Inf, -500, 1, 1, 0, 1, 1, 1),
  m("two-thresh", -500, 500,  0, 0, 0, 1, 1, 1),
  m("two-thresh", 500, Inf,   1, 0, 1, 1, 1, 1))
pred.dt <- data.table(
  problem=unique(model.dt$problem),
  pred.log.lambda=0)

test_that("problem with no thresh is OK", {
  L <- ROChange(model.dt, pred.dt, "problem")
  expect_is(L, "list")
  expect_equal(L$roc$fp, fp <- c(1, 1, 2))
  expect_equal(L$roc$fn, fn <- c(2, 1, 1))
  expect_equal(L$roc$FPR, fp/3)
  expect_equal(L$roc$TPR, 1-fn/3)
})

