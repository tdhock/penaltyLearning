\name{modelSelectionC}
\alias{modelSelectionC}
\title{Exact model selection function}
\description{Given loss.vec L_i, model.complexity K_i, the model selection
function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
the solutions (i, min.lambda, max.lambda) with i being the
solution for every lambda in (min.lambda, max.lambda). This
function uses the linear time algorithm implemented in C code.
This function is mostly meant for internal use -- it is instead
recommended to use modelSelection.}
\usage{modelSelectionC(loss.vec, model.complexity, model.id)}
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
pro <- subset(neuroblastoma$profiles, profile.id==1 & chromosome=="X")
max.segments <- 20
fit <- Segmentor3IsBack::Segmentor(pro$logratio, 2, max.segments)
seg.vec <- 1:max.segments
exact.df <- modelSelectionC(fit@likelihood, seg.vec, seg.vec)
## Solve the optimization using grid search.
L.grid <- with(exact.df,{
  seq(min(max.log.lambda)-1,
      max(min.log.lambda)+1,
      l=100)
})
lambda.grid <- exp(L.grid)
kstar.grid <- sapply(lambda.grid, function(lambda){
  crit <- with(exact.df, model.complexity * lambda + model.loss)
  picked <- which.min(crit)
  exact.df$model.id[picked]
})
grid.df <- data.frame(log.lambda=L.grid, segments=kstar.grid)
library(ggplot2)
## Compare the results.
ggplot()+
  ggtitle("grid search (red) agrees with exact path computation (black)")+
  geom_segment(aes(min.log.lambda, model.id,
                   xend=max.log.lambda, yend=model.id),
               data=exact.df)+
  geom_point(aes(log.lambda, segments),
             data=grid.df, color="red", pch=1)+
  ylab("optimal model complexity (segments)")+
  xlab("log(lambda)")

}
