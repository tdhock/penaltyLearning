\name{modelSelection}
\alias{modelSelection}
\title{Exact model selection function}
\description{Given a data.frame with "loss" column L_i, "model.complexity"
column K_i, model selection function i*(lambda) = argmin_i L_i +
lambda*K_i, compute all of the solutions (i, min.lambda,
max.lambda) with i being the solution for every lambda in
(min.lambda, max.lambda).}
\usage{modelSelection()}

\value{data.frame with a row for each model that can be selected for at
least one lambda value, and the following columns. (min.lambda,
max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
optimal penalty constants, on the original and log scale;
model.complexity are the K_i values; model.id are the model
identifiers (also used for row names); and model.loss are the C_i
values.}

\author{Toby Dylan Hocking}




\examples{

if(require(neuroblastoma) && require(cghseg)){

  data(neuroblastoma, envir=environment())
  pro <- subset(neuroblastoma$profiles, profile.id==1 & chromosome=="X")
  max.segments <- 20
  fit <- cghseg:::segmeanCO(pro$logratio, Kmax=max.segments)
  seg.vec <- 1:max.segments
  exact.df <- modelSelection(fit$J.est, seg.vec, seg.vec)

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

  if(require(ggplot2)){
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
  
}

}
