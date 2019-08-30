\name{featureVector}
\alias{featureVector}
\title{featureVector}
\description{Compute a feature vector of constant length which can be used as
an input for supervised penalty learning. The output is a target
interval of log(penalty) values that achieve minimum incorrect
labels (see \code{\link{targetIntervals}}).}
\usage{featureVector(data.vec)}
\arguments{
  \item{data.vec}{numeric vector of ordered data.}
}

\value{Numeric vector of features.}

\author{Toby Dylan Hocking}




\examples{

data(neuroblastoma, package="neuroblastoma", envir=environment())

one <- subset(neuroblastoma$profiles, profile.id=="1" & chromosome=="1")
(f.vec <- featureVector(one$logratio))

}
