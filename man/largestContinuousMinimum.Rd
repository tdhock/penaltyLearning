\name{largestContinuousMinimum}
\alias{largestContinuousMinimum}
\title{largestContinuousMinimum}
\description{Find the run of minimum cost with the largest size.}
\usage{largestContinuousMinimum(cost, size)}
\arguments{
  \item{cost}{
}
  \item{size}{
}
}



\author{Toby Dylan Hocking}




\examples{

if(require(neuroblastoma) && require(cghseg)){

  data(neuroblastoma, envir=environment())
  pid <- 1
  chr <- 1
  pro <- subset(neuroblastoma$profiles, profile.id==pid & chromosome==chr)
  ann <- subset(neuroblastoma$annotations, profile.id==pid & chromosome==chr)
  max.segments <- 20
  fit <- cghseg:::segmeanCO(pro$logratio, Kmax=max.segments)
  seg.vec <- 1:max.segments
  exact.df <- exactModelSelection(fit$J.est, seg.vec, seg.vec)

}

}
