library(testthat)
context("ROChange-no-thresh")
library(penaltyLearning)

m <- function(problem, min.log.lambda, max.log.lambda, errors, fp, fn, labels, possible.fp, possible.fn){
  data.table(problem, min.log.lambda, max.log.lambda, errors, fp, fn, labels, possible.fp, possible.fn)
}
model.dt <- rbind(
  m("no-thresh", -Inf, Inf,   0, 0, 0, 2, 0, 2),
  m("two-thresh", -Inf, -500, 1, 1, 0, 1, 1, 1),
  m("two-thresh", -500, 500,  0, 0, 0, 1, 1, 1),
  m("two-thresh", 500, Inf,   1, 0, 1, 1, 1, 1))
pred.dt <- data.table(
  problem=c("no-thresh", "two-thresh"),
  pred.log.lambda=0)

test_that("problem with no thresh is OK", {
  L <- ROChange(model.dt, pred.dt, "problem")
  expect_is(L, "list")
})
