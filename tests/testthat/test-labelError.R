library(testthat)
context("labelError")
library(penaltyLearning)
ids.str <- paste(c(1, 4, 6, 8, 10, 11))
someProfiles <- function(all.profiles){
  data.table(all.profiles)[profile.id %in% ids.str, ]
}
data(neuroblastoma, package="neuroblastoma")
profiles <- someProfiles(neuroblastoma$profiles)
labels <- someProfiles(neuroblastoma$annotations)
## Plot labels along with noisy data sets.
breakpoint.colors <- c(
  "breakpoint"="#a445ee",
  "normal"="#f6f4bf")
if(interactive()){
  library(ggplot2)
  ggplot()+
    ggtitle("supervised change-point detection = data + labels")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
    geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation),
                  color="grey",
                  data=labels)+
    scale_fill_manual("label", values=breakpoint.colors)+
    geom_point(aes(position/1e6, logratio),
               data=profiles,
               shape=1)+
    scale_x_continuous(
      "position on chromosome (mega bases)",
      breaks=c(100, 200))
}
problem.list <- split(profiles, profiles[, paste(profile.id, chromosome)])
segs.list <- list()
loss.list <- list()
for(problem.i in seq_along(problem.list)){
  problem.name <- names(problem.list)[[problem.i]]
  cat(sprintf(
    "%4d / %4d problems %s\n",
    problem.i, length(problem.list), problem.name))
  pro <- problem.list[[problem.name]]
  meta <- pro[1, .(profile.id, chromosome)]
  max.segments <- min(nrow(pro), 10)
  fit <- Segmentor3IsBack::Segmentor(
    pro$logratio, model=2, Kmax=max.segments)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
    (pro$position[data.before.change]+pro$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro$position))
    seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
    segs.list[[paste(problem.name, n.segments)]] <- data.table(
      meta,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=seg.mean.vec)
  }
  loss.list[[paste(problem.name, n.segments)]] <- data.table(
    meta,
    n.segments=1:max.segments,
    loss=as.numeric(fit@likelihood))
}
loss <- do.call(rbind, loss.list)
segs <- do.call(rbind, segs.list)

selection <- loss[, {
  penaltyLearning::modelSelection(.SD, "loss", "n.segments")
}, by=.(profile.id, chromosome)]
changes <- segs[1 < start, ]
errors <- labelError(
  selection, labels, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  problem.vars=c("profile.id", "chromosome"))
model.counts <- selection[, list(models=.N), by=.(profile.id, chromosome)]
label.counts <- labels[, list(labels=.N), by=.(profile.id, chromosome)]
labeled.model.counts <- model.counts[label.counts, on=list(
  profile.id, chromosome)]
test_that("label error OK when more models than labels", {
  expect_equal(nrow(errors$model.errors), sum(labeled.model.counts$models))
})

test_that("error for missing columns in targetIntervals", {
  expect_error({
    targetIntervals(selection, problem.vars=c("profile.id", "chromosome"))
  }, "models$errors should be the number of incorrect labels", fixed=TRUE)
})

only10 <- changes[n.segments==10]
test_that("error for missing changes", {
  expect_error({
    labelError(
      selection, labels, only10,
      change.var="chromStart",
      label.vars=c("min", "max"),
      problem.vars=c("profile.id", "chromosome"))
  }, "each model should have a different number of changes")
})

targets <- targetIntervals(
  errors$model.errors,
  problem.vars=c("profile.id", "chromosome"))
test_that("errors are reported", {
  expect_true(is.numeric(targets$errors))
})

test_that("label error fails when more labels than models", {
  expect_error({
    labelError(
      selection, neuroblastoma$annotations, changes,
      change.var="chromStart",
      label.vars=c("min", "max"),
      problem.vars=c("profile.id", "chromosome"))
  }, "some labels have no models")
})

## From example(largestContinuousMinimum).
data(neuroblastoma, package="neuroblastoma", envir=environment())
pro4 <- subset(neuroblastoma$profiles, profile.id==4)
ann4 <- subset(neuroblastoma$annotations, profile.id==4)
label <- function(annotation, min, max){
  data.frame(profile.id=4, chromosome="14", min, max, annotation)
}
ann <- rbind(
  ann4,
  label("1change", 70e6, 80e6),
  label("0changes", 20e6, 60e6))
max.segments <- 20
segs4.list <- list()
selection.list <- list()
for(chr in unique(ann$chromosome)){
  pro <- subset(pro4, chromosome==chr)
  fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
  model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
  selection.df <- modelSelection(model.df, complexity="n.segments")
  selection.list[[chr]] <- data.table(chromosome=chr, selection.df)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
    (pro$position[data.before.change]+pro$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro$position))
    segs4.list[[paste(chr, n.segments)]] <- data.table(
      chromosome=chr,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=fit@parameters[n.segments, 1:n.segments])
  }
}
segs4 <- do.call(rbind, segs4.list)
selection <- do.call(rbind, selection.list)
changes <- segs4[1 < start,]
error.list <- labelError(
  selection, ann, changes,
  problem.vars="chromosome", # for all three data sets.
  model.vars="n.segments", # for changes and selection.
  change.var="chromStart", # column of changes with breakpoint position.
  label.vars=c("min", "max")) # limit of labels in ann.
test_that("same number of model.errors as selection", {
  expect_equal(nrow(error.list$model.errors), nrow(selection))
})

## Trivial edge cases.
## 123456789
##  (-0-]
##      (1]
##  | OK
##   | FP
##      | FP
##       | TP
##        | TP
##         | FN
label <- function(annotation, start, end){
  data.frame(prob="five", start, end, annotation)
}
ann.trivial <- rbind(
  label("1change", 6, 8),
  label("0changes", 2, 6))
models <- data.table(
  prob="five",
  complexity=c(-1, -3, -5, -6))
changes <- data.table(
  prob="five",
  pos=c(1, 7, 1, 6),
  complexity=c(-3, -5, -5, -6))
test_that("labelError throws informative errors", {
  expect_error({
    labelError(models, ann.trivial, changes)
  }, "problem.vars should be a character vector of column names present in models, changes, and labels (ID for separate changepoint detection problems)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob")
  }, "label.vars should be a 2-element character vector of labels column names (start and end of labeled region)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=character())
  }, "label.vars should be a 2-element character vector of labels column names (start and end of labeled region)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("foo", "bar"))
  }, "label.vars should be a 2-element character vector of labels column names (start and end of labeled region)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"))
  }, "change.var should be a column name of changes (position of predicted changepoints)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var=c("foo1"))
  }, "change.var should be a column name of changes (position of predicted changepoints)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var=c("pos", "end"))
  }, "change.var should be a column name of changes (position of predicted changepoints)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var="pos")
  }, "model.vars should be a column name of both models and changes (ID for model complexity, typically the number of changepoints or segments)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("start", "end"),
               change.var="pos")
  }, "model.vars should be a column name of both models and changes (ID for model complexity, typically the number of changepoints or segments)", fixed=TRUE)
  expect_error({
    labelError(models, ann.trivial, changes, problem.vars="prob",
               label.vars=c("end", "start"),
               change.var="pos",
               model.vars="complexity")
  }, "label start must be less than end", fixed=TRUE)
})

trivial.list <- labelError(
  models, ann.trivial, changes,
  problem.vars="prob",
  label.vars=c("start", "end"),
  change.var="pos",
  model.vars="complexity")
test_that("1 TP for complexity=-5, 2 errors for -6", {
  trivial.list$model.errors[, {
    expect_equal(complexity, c(-1, -3, -5, -6))
    expect_equal(labels, c(2, 2, 2, 2))
    expect_equal(errors, c(1, 1, 0, 2))
    expect_equal(possible.fp, c(2, 2, 2, 2))
    expect_equal(fp, c(0, 0, 0, 1))
    expect_equal(possible.fn, c(1, 1, 1, 1))
    expect_equal(fn, c(1, 1, 0, 1))
  }]
})

ann.overlap <- rbind(
  label("1change", 5, 8),
  label("0changes", 2, 6))
test_that("error for overlapping labels", {
  expect_error({
    labelError(
      models, ann.overlap, changes,
      problem.vars="prob",
      label.vars=c("start", "end"),
      change.var="pos",
      model.vars="complexity")
  }, "each label end must be <= next label start", fixed=TRUE)
})
