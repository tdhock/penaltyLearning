library(testthat)
context("peaks")
library(penaltyLearning)
library(data.table)
library(PeakSegDP)
library(PeakError)

data(chr11ChIPseq, package="PeakSegDP")
coverage.list <- with(chr11ChIPseq, split(coverage, coverage$sample.id))
label.list <- with(chr11ChIPseq, split(regions, regions$sample.id))

all.errors.list <- list()
for(sample.id in names(coverage.list)){
  sample.coverage <- coverage.list[[sample.id]]
  sample.labels <- label.list[[sample.id]]
  fit <- PeakSegOptimal::PeakSegPDPAchrom(sample.coverage, 5L)
  sample.peaks <- data.table(fit$segments)[status=="peak", ]
  setkey(sample.peaks, peaks)
  selection <- modelSelection(fit$loss, "PoissonLoss", "peaks")
  sample.errors <- data.table(selection)[, {
    peak.val <- peaks
    model.peaks <- sample.peaks[peaks==peak.val, ]
    model.error <- PeakErrorChrom(model.peaks, sample.labels)
    dt <- data.table(model.error)[, list(
      possible.fp=sum(possible.fp),
      fp=sum(fp),
      possible.fn=sum(possible.tp),
      fn=sum(fn),
      labels=.N,
      errors=sum(fp+fn))]
    data.table(.SD, dt)
  }, by=peaks]
  all.errors.list[[sample.id]] <- data.table(sample.id, sample.errors)
}
all.errors <- do.call(rbind, all.errors.list)

target.dt <- targetIntervals(all.errors, "sample.id")
target.mat <- target.dt[, cbind(min.log.lambda, max.log.lambda)]
feature.mat <- t(sapply(coverage.list, function(df){
  c(n.data=nrow(df),
    max.coverage=max(df$count))
}))
fit <- IntervalRegressionUnregularized(
  feature.mat, target.mat, max.iterations=1e3)
pred.dt <- data.table(
  sample.id=rownames(feature.mat),
  pred.log.lambda=as.numeric(fit$predict(feature.mat)))
roc <- ROChange(all.errors, pred.dt, "sample.id")

test_that("perfect prediction train error for peak model", {
  expect_equal(roc$auc, 1)
})
test_that("TPR=1 is achieved but FPR=1 is not", {
  expect_gt(nrow(roc$roc[TPR==1,]), 0)
  expect_equal(nrow(roc$roc[FPR==1,]), 0)
})

## ggplot()+
##   geom_polygon(aes(FPR, TPR),
##                color="black",
##                fill="grey",
##                data=roc$auc.polygon)
