\name{coef.IntervalRegression}
\alias{coef.IntervalRegression}
\title{coef IntervalRegression}
\description{Get the learned coefficients of an IntervalRegression model.}
\usage{\method{coef}{IntervalRegression}(object, 
    ...)}
\arguments{
  \item{object}{
}
  \item{\dots}{
}
}

\value{numeric matrix [features x regularizations] of learned weights (on
the original feature scale), can be used for prediction via
cbind(1,features) \%*\% weights.}

\author{Toby Dylan Hocking}





