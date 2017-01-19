\name{neuroblastomaProcessed}
\alias{neuroblastomaProcessed}
\docType{data}
\title{
  Processed neuroblastoma data set with features and targets
}
\description{
  Features are inputs and targets are outputs for
  penalty learning functions like penaltyLearning::IntervalRegressionCV.
  data(neuroblastoma, package="neuroblastoma") was processed by
  computing optimal Gaussian segmentation models from 1 to 20 segments
  (cghseg:::segmeanCO or Segmentor3IsBack::Segmentor), then label error
  was computed using neuroblastoma$annotations
  (penaltyLearning::labelError), then target intervals were
  computed (penaltyLearning::targetInterval). Features were also
  computed based on neuroblastoma$profiles.
}
\usage{data("neuroblastomaProcessed")}
\format{
List of two matrices: feature.mat is n.observations x n.features, and
  target.mat is n.observations x 2, where n.observations=3418 and
  n.features=117.
}
\keyword{datasets}
