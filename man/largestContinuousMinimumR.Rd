\name{largestContinuousMinimumR}
\alias{largestContinuousMinimumR}
\title{largestContinuousMinimumR}
\description{Find the run of minimum cost with the largest size.}
\usage{largestContinuousMinimumR(cost, size)}
\arguments{
  \item{cost}{numeric vector of cost values.}
  \item{size}{numeric vector of interval size values.}
}

\value{Integer vector length 2 (start and end of target interval relative
to cost and size).}

\author{Toby Dylan Hocking}




\examples{

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
segs.list <- list()
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
    segs.list[[paste(chr, n.segments)]] <- data.table(
      chromosome=chr,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=fit@parameters[n.segments, 1:n.segments])
  }
}
segs <- do.call(rbind, segs.list)
selection <- do.call(rbind, selection.list)

changes <- segs[1 < start,]
error.list <- labelError(
  selection, ann, changes,
  problem.vars="chromosome", # for all three data sets.
  model.vars="n.segments", # for changes and selection.
  change.var="chromStart", # column of changes with breakpoint position.
  label.vars=c("min", "max")) # limit of labels in ann.

one.problem.error <- error.list$model.errors[chromosome=="14", ]
indices <- one.problem.error[, largestContinuousMinimumR(
  errors, max.log.lambda-min.log.lambda)]
one.problem.error[indices[["start"]]:indices[["end"]],]

}
