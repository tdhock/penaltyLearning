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
setkey(model.counts, profile.id, chromosome)
setkey(label.counts, profile.id, chromosome)
labeled.model.counts <- model.counts[label.counts]
test_that("label error OK when more models than labels", {
  expect_equal(nrow(errors$model.errors), sum(labeled.model.counts$models))
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
