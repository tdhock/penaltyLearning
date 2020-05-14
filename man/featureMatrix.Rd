\name{featureMatrix}
\alias{featureMatrix}
\title{featureMatrix}
\description{Compute a feature matrix (segmentation problems x features).}
\usage{featureMatrix(data.sequences, 
    problem.vars, data.var)}
\arguments{
  \item{data.sequences}{data.frame of sorted sequences of data to segment.}
  \item{problem.vars}{character vector of columns of \code{data.sequences} to treat as
segmentation problem IDs.}
  \item{data.var}{character vector of length 1 (column of \code{data.sequences} to treat as
data to segment).}
}

\value{Numeric feature matrix. Some entries may be missing or infinite;
these columns should be removed before model training.}

\author{Toby Dylan Hocking}




\examples{

test.df <- data.frame(
  id=rep(1:2, each=10),
  x=rnorm(20))
penaltyLearning::featureMatrix(test.df, "id", "x")
if(requireNamespace("neuroblastoma")){
  data(neuroblastoma, package="neuroblastoma", envir=environment())
  one <- subset(neuroblastoma$profiles, profile.id \%in\% c(1,2))
  f.mat <- penaltyLearning::featureMatrix(
    one, c("profile.id", "chromosome"), "logratio")
}

}
