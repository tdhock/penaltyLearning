\name{targetIntervals}
\alias{targetIntervals}
\title{Compute target intervals}
\description{Compute target intervals of log(penalty) values that result in
predicted changepoint models with minimum incorrect labels.
Use this function after labelError, and before IntervalRegression*.}
\usage{targetIntervals(models, problem.vars)}
\arguments{
  \item{models}{data.table with columns errors, min.log.lambda, max.log.lambda,
typically labelError()$model.errors.}
  \item{problem.vars}{character: column names used to identify data set / segmentation
problem.}
}

\value{data.table with columns problem.vars, one row for each
segmentation problem. The "min.log.lambda", and "max.log.lambda"
columns give the largest interval of log(penalty) values which
results in the minimum incorrect labels for that problem. This can
be used to create the target.mat parameter of the
IntervalRegression* functions.}

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
targetIntervals(error.list$model.errors, "chromosome")

}
