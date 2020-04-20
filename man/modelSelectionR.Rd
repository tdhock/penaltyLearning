\name{modelSelectionR}
\alias{modelSelectionR}
\title{Exact model selection function}
\description{Given \code{loss.vec} L_i, \code{model.complexity} K_i, the model selection
function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
the solutions (i, min.lambda, max.lambda) with i being the
solution for every lambda in (min.lambda, max.lambda). This
function uses the quadratic time algorithm implemented in R code.
This function is mostly meant for internal use -- it is instead
recommended to use \code{\link{modelSelection}}.}
\usage{modelSelectionR(loss.vec, 
    model.complexity, 
    model.id)}
\arguments{
  \item{loss.vec}{numeric vector: loss L_i}
  \item{model.complexity}{numeric vector: model complexity K_i}
  \item{model.id}{vector: indices i}
}

\value{data.frame with a row for each model that can be selected for at
least one lambda value, and the following columns. (min.lambda,
max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
optimal penalty constants, on the original and log scale;
\code{model.complexity} are the K_i values; \code{model.id} are the model
identifiers (also used for row names); and model.loss are the C_i
values.}

\author{Toby Dylan Hocking}




\examples{

if(interactive()){
  library(penaltyLearning)
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  one <- subset(neuroblastoma$profiles, profile.id==599 & chromosome=="14")
  max.segments <- 1000
  fit <- Segmentor3IsBack::Segmentor(one$logratio, model=2, Kmax=max.segments)
  lik.df <- data.frame(lik=fit@likelihood, segments=1:max.segments)
  times.list <- list()
  for(n.segments in seq(10, max.segments, by=10)){
    some.lik <- lik.df[1:n.segments,]
    some.times <- microbenchmark::microbenchmark(
      R=pathR <- with(some.lik, modelSelectionR(lik, segments, segments)),
      C=pathC <- with(some.lik, modelSelectionC(lik, segments, segments)),
      times=5)
    times.list[[paste(n.segments)]] <- data.frame(n.segments, some.times)
  }
  times <- do.call(rbind, times.list)
  ## modelSelectionR and modelSelectionC should give identical results.
  identical(pathR, pathC)
  ## However, modelSelectionC is much faster (linear time complexity)
  ## than modelSelectionR (quadratic time complexity).
  library(ggplot2)
  library(data.table)
  times.dt <- data.table(times)
  times.dt[, seconds := time/1e9]
  gg <- ggplot()+
    geom_point(aes(
      n.segments, seconds, color=expr),
      data=times.dt,
      shape=1)
  gg+
    scale_x_log10()+
    scale_y_log10()
  R.times <- times.dt[expr=="R"]
  R.times[, log.seconds := log(seconds)]
  R.times[, log.segs := log(n.segments)]
  R.times[, segs.squared := n.segments*n.segments]
  thresh <- 700
  big.segs <- R.times[thresh < n.segments]
  (log.fit <- lm(log.seconds ~ log.segs, data=big.segs))
  (quad.fit <- lm(seconds ~ n.segments + segs.squared, data=R.times))
  big.segs[, pred.seconds := exp(predict(log.fit))]
  R.times[, pred.seconds := predict(quad.fit)]
  pred.dt <- rbind(
    big.segs[, data.table(model="log", expr="R", pred.seconds, n.segments)],
    R.times[, data.table(model="quad", expr="R", pred.seconds, n.segments)])
  ggfit <- gg+
    geom_line(aes(
      n.segments, pred.seconds, linetype=model),
      data=pred.dt)
  ggfit+
    scale_x_log10()+
    scale_y_log10()
  extrapolate.dt <- data.table(
    n.segments=c(100000, 250000, 500000, 1e6))
  extrapolate.dt[, log.segs := log(n.segments)]
  extrapolate.dt[, segs.squared := n.segments*n.segments]
  pred.seconds.mat <- rbind(
    quad=predict(quad.fit, extrapolate.dt),
    log=exp(predict(log.fit, extrapolate.dt)))
  colnames(pred.seconds.mat) <- extrapolate.dt$n.segments
  (pred.hours.mat <- pred.seconds.mat/60/60)
}

}
