library(testthat)
context("modelSelection")
library(penaltyLearning)
data(oneSkip)

test_that("no error for odd model/loss values", {
  loss.df <- data.frame(
    complexity=0:10,
    loss=c(10:6, 0, 5:1))
  selection.df <- modelSelection(loss.df)
  expect_equal(selection.df$complexity, c(5, 0))
})

test_that("output intervals computed correctly", {
  df <- with(oneSkip$input, modelSelectionC(error, segments, peaks))
  expect_identical(df$model.complexity, oneSkip$output$model.complexity)
  expect_identical(df$min.lambda, oneSkip$output$min.lambda)
})

unsorted <- rbind(
  data.frame(segments=c(3, 13), peaks=c(1, 6), error=-3e6),
  oneSkip$input[c(6,4,5,5,3,1,2),])
test_that("modelSelectionC errors for unsorted data", {
  expect_error({
    with(unsorted, modelSelectionC(error, segments, peaks))
  })
})
test_that("modelSelection works for unsorted data", {
  df <- modelSelection(unsorted, "error", "segments")
  expect_identical(df$segments, oneSkip$output$model.complexity)
  expect_identical(df$min.lambda, oneSkip$output$min.lambda)
})

test_that("error when models is not DF",{
  expect_error({
    modelSelection(1L)
  }, "models must be data.frame with at least one row and numeric columns models[[complexity]] and models[[loss]]", fixed=TRUE)
})

test_that("error when models has 0 rows",{
  expect_error({
    modelSelection(data.frame())
  }, "models must be data.frame with at least one row and numeric columns models[[complexity]] and models[[loss]]", fixed=TRUE)
})

test_that("error when models has 1 missing character row",{
  expect_error({
    modelSelection(data.frame(loss=NA_character_, complexity=NA_character_))
  }, "models must be data.frame with at least one row and numeric columns models[[complexity]] and models[[loss]]", fixed=TRUE)
})

test_that("error when models has 1 missing loss",{
  expect_error({
    modelSelection(data.frame(loss=NA_real_, complexity=1))
  }, "which are not missing/NA", fixed=TRUE)
})

test_that("error when models has 1 missing complexity",{
  expect_error({
    modelSelection(data.frame(loss=1, complexity=NA_real_))
  }, "which are not missing/NA", fixed=TRUE)
})

test_that("one model is fine",{
  one <- modelSelection(data.frame(loss=1, complexity=1, foo="bar"))
  expect_identical(paste(one$foo), "bar")
  expect_identical(one$min.lambda, 0)
  expect_identical(one$max.lambda, Inf)
  expect_identical(one$min.log.lambda, -Inf)
  expect_identical(one$max.log.lambda, Inf)
  expect_identical(one$loss, 1)
  expect_identical(one$complexity, 1)
})

test_that("error for bad column names", {
  expect_error({
    modelSelection(loss=NULL)
  }, "loss must be a column name of models")
  expect_error({
    modelSelection(loss=c())
  }, "loss must be a column name of models")
  expect_error({
    modelSelection(loss=c("foo", "bar"))
  }, "loss must be a column name of models")
})

## trivial.
loss.vec <- c(5,4,4)
model.complexity <- c(1,2,3)
n.models <- 3
test_that("loss not decreasing error in C code", {
  expect_error({
    .C(
      "modelSelection_interface",
      loss.vec=as.double(loss.vec),
      model.complexity=as.double(model.complexity),
      n.models=as.integer(n.models),
      after.vec=integer(n.models),
      lambda.vec=double(n.models),
      PACKAGE="penaltyLearning")
  }, "loss not decreasing")
})
test_that("loss not decreasing error in C code Fwd", {
  expect_error({
    .C(
      "modelSelectionFwd_interface",
      loss.vec=as.double(loss.vec),
      model.complexity=as.double(model.complexity),
      n.models=as.integer(n.models),
      after.vec=integer(n.models),
      lambda.vec=double(n.models),
      iterations=integer(n.models),
      PACKAGE="penaltyLearning")
  }, "loss not decreasing")
})

loss.vec <- c(5,4,3)
model.complexity <- c(1,2,2)
n.models <- 3
test_that("complexity not increasing error in C code", {
  expect_error({
    .C(
      "modelSelection_interface",
      loss.vec=as.double(loss.vec),
      model.complexity=as.double(model.complexity),
      n.models=as.integer(n.models),
      after.vec=integer(n.models),
      lambda.vec=double(n.models),
      PACKAGE="penaltyLearning")
  }, "complexity not increasing")
})
test_that("complexity not increasing error in C code Fwd", {
  expect_error({
    .C(
      "modelSelectionFwd_interface",
      loss.vec=as.double(loss.vec),
      model.complexity=as.double(model.complexity),
      n.models=as.integer(n.models),
      after.vec=integer(n.models),
      lambda.vec=double(n.models),
      iterations=integer(n.models),
      PACKAGE="penaltyLearning")
  }, "complexity not increasing")
})

## synthetic data from paper.
N <- 5
t <- 1:N
test_that("2N-3 iterations for worst case synthetic loss values", {
  df <- data.frame(loss=N-t+(1 < t & t < N)/2, complexity=t)
  result <- .C(
    "modelSelectionFwd_interface",
    loss=as.double(df$loss),
    complexity=as.double(df$complexity),
    N=as.integer(nrow(df)),
    models=integer(nrow(df)),
    breaks=double(nrow(df)),
    evals=integer(nrow(df)),
    PACKAGE="penaltyLearning")
  expect_equal(result$evals, c(0, 1, 2, 2, 2))
})

test_that("2N-3 iterations for another worst case", {
  df <- data.frame(loss=N-t+N*(t!=N), complexity=t)
  result <- .C(
    "modelSelectionFwd_interface",
    loss=as.double(df$loss),
    complexity=as.double(df$complexity),
    N=as.integer(nrow(df)),
    models=integer(nrow(df)),
    breaks=double(nrow(df)),
    evals=integer(nrow(df)),
    PACKAGE="penaltyLearning")
  expect_equal(result$evals, c(0, 1, 2, 2, 2))
})

test_that("2N-3 iterations for simple worst case", {
  df <- data.frame(loss=N-t, complexity=t)
  result <- .C(
    "modelSelectionFwd_interface",
    loss=as.double(df$loss),
    complexity=as.double(df$complexity),
    N=as.integer(nrow(df)),
    models=integer(nrow(df)),
    breaks=double(nrow(df)),
    evals=integer(nrow(df)),
    PACKAGE="penaltyLearning")
  expect_equal(result$evals, c(0, 1, 2, 2, 2))
})

test_that("N-1 iterations for best case", {
  df <- data.frame(loss=N-log(t), complexity=t)
  if(FALSE){
    library(ggplot2)
    ggplot()+
      geom_abline(aes(
        slope=complexity, intercept=loss),
        data=df)+
      xlim(0, 10)+
      ylim(0, 10)
  }
  result <- .C(
    "modelSelectionFwd_interface",
    loss=as.double(df$loss),
    complexity=as.double(df$complexity),
    N=as.integer(nrow(df)),
    models=integer(nrow(df)),
    breaks=double(nrow(df)),
    evals=integer(nrow(df)),
    PACKAGE="penaltyLearning")
  expect_equal(result$evals, c(0, 1, 1, 1, 1))
})

test_that("N-1 iterations for another best case", {
  df <- data.frame(loss=N-round(sqrt(t), 2), complexity=t)
  result <- .C(
    "modelSelectionFwd_interface",
    loss=as.double(df$loss),
    complexity=as.double(df$complexity),
    N=as.integer(nrow(df)),
    models=integer(nrow(df)),
    breaks=double(nrow(df)),
    evals=integer(nrow(df)),
    PACKAGE="penaltyLearning")
  if(FALSE){
    library(ggplot2)
    ggplot()+
      geom_abline(aes(
        slope=complexity, intercept=loss),
        data=df)+
      xlim(0, 1)+
      ylim(3, 5)+
      geom_point(aes(
        penalty, cost),
        data=with(result, data.frame(
          penalty=breaks, cost=loss+complexity*breaks)))
  }
  expect_equal(result$evals, c(0, 1, 1, 1, 1))
})
