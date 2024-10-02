\name{labelError}
\alias{labelError}
\title{Compute incorrect labels}
\description{Compute incorrect \code{labels} for several change-point detection
problems and \code{models}. Use this function after having computed
changepoints, loss values, and model selection functions
(see \code{\link{modelSelection}}). The next step after labelError is typically
computing target intervals of log(penalty) values that predict
changepoints with minimum incorrect \code{labels} for each problem (see
\code{\link{targetIntervals}}).}
\usage{labelError(models, labels, 
    changes, change.var = "chromStart", 
    label.vars = c("min", 
        "max"), model.vars = "n.segments", 
    problem.vars = character(0), 
    annotations = change.labels)}
\arguments{
  \item{models}{data.frame with one row per (problem,model) combination, typically
the output of modelSelection(...). There is a row for each
changepoint model that could be selected for a particular
segmentation problem. There should be columns \code{problem.vars} (for
problem ID) and \code{model.vars} (for model complexity).}
  \item{labels}{data.frame with one row per (problem,region). Each label defines a
region in a particular segmentation problem, and a range of
predicted changepoints which are consistent in that region. There
should be a column "annotation" with takes one of the
corresponding values in the annotation column of \code{\link{change.labels}}
(used to determine the range of predicted changepoints which are
consistent). There should also be a columns \code{problem.vars} (for
problem ID) and \code{label.vars} (for region start/end).}
  \item{changes}{data.frame with one row per (problem,model,change), for each
predicted changepoint (in each model and segmentation
problem). Should have columns \code{problem.vars} (for problem ID),
\code{model.vars} (for model complexity), and \code{change.var} (for changepoint
position).}
  \item{change.var}{character(length=1): column name of predicted change-point
position in \code{labels}. The default "chromStart" is useful for genomic
data with segment start/end positions stored in columns named
chromStart/chromEnd. A predicted changepoint at position X is
interpreted to mean a changepoint between X and X+1.}
  \item{label.vars}{character(length=2): column names of start and end positions of
\code{labels}, in same units as change-point positions. The default is
c("min", "max"). Labeled regions are (start,end] -- open on the
left and closed on the right, so for example a 0changes annotation
between start=10 and end=20 means that any predicted changepoint
at 11, ..., 20 is a false positive.}
  \item{model.vars}{character: column names used to identify model complexity. The
default "n.segments" is for change-point \code{models} such as in the
jointseg and changepoint packages.}
  \item{problem.vars}{character: column names used to identify data set / segmentation
problem, should be present in all three data tables (\code{models},
\code{labels}, \code{changes}).}
  \item{annotations}{data.table with columns annotation, min.changes, max.changes,
possible.fn, possible.fp which is joined to \code{labels} in order to
determine how to compute false positives and false negatives for
each annotation.}
}

\value{list of two data.tables: label.errors has one row for every
combination of \code{models} and \code{labels}, with status column that
indicates whether or not that model commits an error in that
particular label; model.errors has one row per model, with columns
for computing target intervals and ROC curves (see \code{\link{targetIntervals}}
and \code{\link{ROChange}}).}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

label <- function(annotation, min, max){
  data.frame(profile.id=4, chrom="chr14", min, max, annotation)
}
label.df <- rbind(
  label("1change", 70e6, 80e6),
  label("0changes", 20e6, 60e6))
model.df <- data.frame(chrom="chr14", n.segments=1:3)
change.df <- data.frame(chrom="chr14", rbind(
  data.frame(n.segments=2, changepoint=75e6),
  data.frame(n.segments=3, changepoint=c(75e6, 50e6))))
penaltyLearning::labelError(
  model.df, label.df, change.df,
  problem.vars="chrom", # for all three data sets.
  model.vars="n.segments", # for changes and selection.
  change.var="changepoint", # column of changes with breakpoint position.
  label.vars=c("min", "max")) # limit of labels in ann.

}
