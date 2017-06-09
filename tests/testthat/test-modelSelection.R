library(testthat)
context("modelSelection")
library(penaltyLearning)
data(oneSkip)

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

## This is a data set with n=9 data points to segments, but only 8
## unique data points.
small <- subset(neuroblastoma$profiles, profile.id==203 & chromosome=="Y")
max.segments <- 9
fit <- Segmentor(small$logratio, model=2, Kmax=max.segments)
rss.vec <- rep(NA, max.segments)
for(n.segments in 1:max.segments){
  end <- fit@breaks[n.segments, 1:n.segments]
  seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
  data.before.change <- end[-n.segments]
  data.after.change <- data.before.change+1
  pos.before.change <- as.integer(
    (small$position[data.before.change]+small$position[data.after.change])/2)
  start <- c(1, data.after.change)
  data.mean.vec <- rep(seg.mean.vec, end-start+1)
  rss.vec[n.segments] <- sum((small$logratio-data.mean.vec)^2)
}
loss.df <- data.frame(
  n.segments=1:max.segments,
  lik=fit@likelihood,
  loss=rss.vec)
test_that("biggest n.segments=8 for rss loss", {
  selection.df <- modelSelection(
    loss.df, complexity="n.segments")
  expect_equal(max(selection.df$n.segments), 8)
})
test_that("biggest n.segments=8 for lik loss", {
  selection.df <- modelSelection(
    loss.df, complexity="n.segments", loss="lik")
  expect_equal(max(selection.df$n.segments), 8)
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
