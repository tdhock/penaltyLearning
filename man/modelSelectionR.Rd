\name{modelSelectionR}
\alias{modelSelectionR}
\title{Exact model selection function}
\description{Given loss.vec L_i, model.complexity K_i, the model selection
function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
the solutions (i, min.lambda, max.lambda) with i being the
solution for every lambda in (min.lambda, max.lambda). This
function uses the quadratic time algorithm implemented in R code.
This function is mostly meant for internal use -- it is instead
recommended to use modelSelection.}
\usage{modelSelectionR(loss.vec, model.complexity, model.id)}
\arguments{
  \item{loss.vec}{numeric vector: loss L_i}
  \item{model.complexity}{numeric vector: model complexity K_i}
  \item{model.id}{vector: indices i}
}

\value{data.frame with a row for each model that can be selected for at
least one lambda value, and the following columns. (min.lambda,
max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
optimal penalty constants, on the original and log scale;
model.complexity are the K_i values; model.id are the model
identifiers (also used for row names); and model.loss are the C_i
values.}

\author{Toby Dylan Hocking}




\examples{

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
ggplot()+
  geom_point(aes(n.segments, time/1e9, color=expr), data=times)

}
