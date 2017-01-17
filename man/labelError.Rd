\name{labelError}
\alias{labelError}
\title{Compute incorrect labels}
\description{Compute incorrect labels for several change-point detection
problems and models.}
\usage{labelError(models, labels, changes, change.var = "chromStart", 
    label.vars = c("min", "max"), model.vars = "n.segments", 
    problem.vars = character(0))}
\arguments{
  \item{models}{data.frame with one row per (problem,model) combination.}
  \item{labels}{data.frame with one row per (problem,label).}
  \item{changes}{data.frame with one row per (problem,model,change).}
  \item{change.var}{character(length=1): column name of predicted change-point
position (refers to the changes argument). The default
"chromStart" is useful for genomic data with segment start/end
positions stored in columns named chromStart/chromEnd.}
  \item{label.vars}{character(length=2): column names of start and end positions of
labeled regions, in same units as change-point positions (refers
to the labels argument). The default is c("min", "max").}
  \item{model.vars}{character: column names used to identify model complexity. The
default "n.segments" is for change-point models such as in the
Segmentor3IsBack and cghseg packages.}
  \item{problem.vars}{character: column names used to identify data set / segmentation
problem. }
}

\value{list of two data.tables: label.errors has one row for every
combination of models and labels, with status column that
indicates whether or not that model commits an error in that
particular label; model.errors has one row per row of models, with
columns for computing error and ROC curves.}

\author{Toby Dylan Hocking}




\examples{

library(penaltyLearning)
data(neuroblastoma, package="neuroblastoma", envir=environment())
pro4 <- subset(neuroblastoma$profiles, profile.id==4)
ann4 <- subset(neuroblastoma$annotations, profile.id==4)
label <- function(annotation, min, max){
  data.table(profile.id=4, chromosome="14", min, max, annotation)
}
ann <- rbind(
    ann4,
    label("1change", 70e6, 80e6),
    label("0changes", 20e6, 60e6))
max.segments <- 5
segs.list <- list()
models.list <- list()
for(chr in unique(ann$chromosome)){
  pro <- subset(pro4, chromosome==chr)
  fit <- Segmentor3IsBack::Segmentor(pro$logratio, model=2, Kmax=max.segments)
  model.df <- data.frame(loss=fit@likelihood, n.segments=1:max.segments)
  models.list[[chr]] <- data.table(chromosome=chr, model.df)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
    (pro$position[data.before.change]+pro$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro$position))
    segs.list[[paste(chr, n.segments)]] <- data.table(
        chromosome=chr,
        n.segments,
        start,
        end,
        chromStart,
        chromEnd,
        mean=fit@parameters[n.segments, 1:n.segments])
  }
}
segs <- do.call(rbind, segs.list)
models <- do.call(rbind, models.list)

changes <- segs[1 < start,]
error.list <- labelError(
    models, ann, changes,
    problem.vars="chromosome", # for all three data sets.
    model.vars="n.segments", # for changes and selection.
    change.var="chromStart", # column of changes with breakpoint position.
    label.vars=c("min", "max")) # limit of labels in ann.

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(n.segments ~ chromosome, scales="free", space="free")+
  scale_x_continuous(breaks=c(100, 200))+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                                 "false negative"=3,
                                 "false positive"=1))+
  scale_fill_manual("label", values=change.colors)+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6),
                color="grey",
                fill=NA,
                data=error.list$label.errors)+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6,
                    fill=annotation, linetype=status),
                data=error.list$label.errors)+
  geom_point(aes(position/1e6, logratio),
             data=subset(pro4, chromosome \%in\% ann$chromosome),
             shape=1)+
  geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
               data=segs,
               color="green",
               size=1)+
  geom_vline(aes(xintercept=chromStart/1e6),
             data=changes,
             linetype="dashed",
             color="green")

}
