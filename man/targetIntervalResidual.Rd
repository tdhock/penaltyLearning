\name{targetIntervalResidual}
\alias{targetIntervalResidual}
\title{targetIntervalResidual}
\description{Compute residual of predicted penalties with respect to target
intervals. This function is useful for visualizing the errors in a
plot of log(penalty) versus a feature.}
\usage{targetIntervalResidual(target.mat, pred)}
\arguments{
  \item{target.mat}{n x 2 numeric matrix: target intervals of log(penalty) values that
yield minimal incorrect labels.}
  \item{pred}{numeric vector: predicted log(penalty) values.}
}

\value{numeric vector of n residuals. Predictions that are too high
(above target.mat[,2]) get positive residuals (too few
changepoints), and predictions that are too low (below
target.mat[,1]) get negative residuals.}

\author{Toby Dylan Hocking}




\examples{

library(penaltyLearning)
data(neuroblastoma, package="neuroblastoma", envir=environment())
pid <- c(4, 520)
pro <- subset(neuroblastoma$profiles, profile.id \%in\% pid)
ann <- subset(neuroblastoma$annotations, profile.id \%in\% pid)
max.segments <- 20
segs.list <- list()
selection.list <- list()
for(i in 1:nrow(ann)){
  a <- ann[i,]
  id.str <- paste(a$chromosome, a$profile.id)
  pro.chr <- subset(pro, chromosome==a$chromosome & profile.id==a$profile.id)
  fit <- Segmentor3IsBack::Segmentor(
    pro.chr$logratio, model=2, Kmax=max.segments)
  model.dt <- data.table(
    lik=as.numeric(fit@likelihood),
    n.segments=1:max.segments,
    loss=NA_real_)
  for(n.segments in 1:max.segments){
    end <- fit@breaks[n.segments, 1:n.segments]
    data.before.change <- end[-n.segments]
    data.after.change <- data.before.change+1
    pos.before.change <- as.integer(
    (pro.chr$position[data.before.change]+
     pro.chr$position[data.after.change])/2)
    start <- c(1, data.after.change)
    chromStart <- c(pro.chr$position[1], pos.before.change)
    chromEnd <- c(pos.before.change, max(pro.chr$position))
    seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
    data.mean.vec <- rep(seg.mean.vec, end-start+1)
    residual.vec <- pro.chr$logratio-data.mean.vec
    model.dt[n.segments, loss := sum(residual.vec * residual.vec)]
    segs.list[[paste(id.str, n.segments)]] <- data.table(
      profile.id=a$profile.id,
      chromosome=a$chromosome,
      n.segments,
      start,
      end,
      chromStart,
      chromEnd,
      mean=seg.mean.vec)
  }
  selection.dt <- modelSelection(model.dt, complexity="n.segments")
  selection.list[[id.str]] <- data.table(
    profile.id=a$profile.id,
    chromosome=a$chromosome,
    selection.dt)
}
segs <- do.call(rbind, segs.list)
selection <- do.call(rbind, selection.list)
changes <- segs[1 < start,]
error.list <- labelError(
  selection, ann, changes,
  problem.vars=c("profile.id", "chromosome"),
  model.vars="n.segments", 
  change.var="chromStart", 
  label.vars=c("min", "max")) 
target.dt <- targetIntervals(
  error.list$model.errors, c("profile.id", "chromosome"))
## The BIC model selection criterion is lambda = log(n), where n is
## the number of data points to segment. This implies log(lambda) =
## log(log(n)), which is a feature we can compute:
feature.dt <- data.table(pro)[, list(
  loglog.n=log(log(.N))
), by=list(profile.id, chromosome)]
pred.dt <- feature.dt[target.dt, on=list(profile.id, chromosome)]
pred.dt[, pred.log.lambda := loglog.n ]
pred.dt[, residual := targetIntervalResidual(
  cbind(min.log.lambda, max.log.lambda),
  pred.log.lambda)]

library(ggplot2)
limits.dt <- pred.dt[, data.table(
  loglog.n,
  log.penalty=c(min.log.lambda, max.log.lambda),
  limit=rep(c("min", "max"), each=.N))][is.finite(log.penalty)]
ggplot()+
  geom_abline(slope=1, intercept=0)+
  geom_point(aes(
    loglog.n,
    log.penalty,
    fill=limit),
    data=limits.dt,
    shape=21)+
  geom_segment(aes(
    loglog.n, pred.log.lambda,
    xend=loglog.n, yend=pred.log.lambda-residual),
    data=pred.dt,
    color="red")+
  scale_fill_manual(values=c(min="white", max="black"))

}
