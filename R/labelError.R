labelError <- structure(function # Compute incorrect labels
### Compute incorrect labels for several change-point detection
### problems and models.
(models,
### data.frame with one row per (problem,model) combination.
  labels,
### data.frame with one row per (problem,label).
  changes,
### data.frame with one row per (problem,model,change).
  change.var="chromStart",
### character(length=1): column name of predicted change-point
### position (refers to the changes argument). The default
### "chromStart" is useful for genomic data with segment start/end
### positions stored in columns named chromStart/chromEnd.
  label.vars=c("min", "max"),
### character(length=2): column names of start and end positions of
### labeled regions, in same units as change-point positions (refers
### to the labels argument). The default is c("min", "max").
  model.vars="n.segments",
### character: column names used to identify model complexity. The
### default "n.segments" is for change-point models such as in the
### Segmentor3IsBack and cghseg packages.
  problem.vars=character(0)
### character: column names used to identify data set / segmentation
### problem. 
){
  stopifnot(is.character(problem.vars))
  stopifnot(is.character(model.vars))
  stopifnot(is.character(change.var))
  stopifnot(is.character(label.vars))
  stopifnot(length(change.var)==1)
  stopifnot(length(label.vars)==2)
  stopifnot(is.data.frame(models))
  stopifnot(is.data.frame(labels))
  stopifnot(is.data.frame(changes))
  if(length(problem.vars)==0){
    stop("Need at least one column name in problem.vars")
  }
  if(length(model.vars)==0){
    stop("Need at least one column name in model.vars")
  }
  stopifnot(label.vars %in% names(labels))
  stopifnot(change.var %in% names(changes))
  stopifnot(problem.vars %in% names(changes))
  stopifnot(problem.vars %in% names(labels))
  stopifnot(problem.vars %in% names(models))
  stopifnot(model.vars %in% names(models))
  stopifnot(model.vars %in% names(changes))
  new.key <- paste0(change.var, ".after")
  if(new.key %in% names(changes)){
    stop("changes should not have a column named ", new.key)
  }
  labels.dt <- data.table(labels)
  if("weight" %in% names(labels.dt)){
    stopifnot(is.numeric(labels.dt$weight))
    stopifnot(0 < labels.dt$weight)
  }else{
    labels.dt$weight <- 1
  }
  setkey(labels.dt, annotation)
  labels.info <- change.labels[labels.dt]
  stopifnot(nrow(labels.info)==nrow(labels.dt))
  setkeyv(labels.info, problem.vars)
  models.dt <- data.table(models)
  setkeyv(models.dt, problem.vars)
  model.labels <- labels.info[models.dt, allow.cartesian=TRUE]
  changes.dt <- data.table(changes)
  changes.dt[[new.key]] <- changes.dt[[change.var]]+1
  changes.key <- c(problem.vars, model.vars, change.var, new.key)
  setkeyv(changes.dt, changes.key)
  labels.key <- c(problem.vars, model.vars, label.vars)
  setkeyv(model.labels, labels.key)
  over.dt <- foverlaps(changes.dt, model.labels, nomatch=0L)
  long.key <- c(
    labels.key, "annotation",
    "min.changes", "max.changes", "possible.fp", "possible.fn", "weight")
  setkeyv(over.dt, long.key)
  setkeyv(model.labels, long.key)
  changes.per.label <- over.dt[model.labels, list(
    pred.changes=.N
    ), by=.EACHI]
  changes.per.label[, fp := ifelse(max.changes < pred.changes, weight, 0)]
  changes.per.label[, fn := ifelse(pred.changes < min.changes, weight, 0)]
  changes.per.label[, status := ifelse(
    fp, "false positive", ifelse(
      fn, "false negative", "correct"))]
  setkeyv(models.dt, c(problem.vars, model.vars))
  setkeyv(changes.per.label, c(problem.vars, model.vars))
  error.totals <- changes.per.label[models.dt, list(
    possible.fp=sum(possible.fp*weight),
    fp=sum(fp),
    possible.fn=sum(possible.fn*weight),
    fn=sum(fn),
    labels=sum(weight),
    errors=sum(fp+fn)),
    by=.EACHI][models.dt]
  list(model.errors=error.totals, label.errors=changes.per.label)
### list of two data.tables: label.errors has one row for every
### combination of models and labels, with status column that
### indicates whether or not that model commits an error in that
### particular label; model.errors has one row per row of models, with
### columns for computing error and ROC curves.
}, ex=function(){
  
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
               data=subset(pro4, chromosome %in% ann$chromosome),
               shape=1)+
    geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
                 data=segs,
                 color="green",
                 size=1)+
    geom_vline(aes(xintercept=chromStart/1e6),
               data=changes,
               linetype="dashed",
               color="green")
  
})

### Describe an annotated region label for supervised change-point detection.
changeLabel <- function(annotation, min.changes, max.changes, color){
  data.table(annotation, min.changes, max.changes, color)
}
change.labels <- rbind(
  changeLabel("breakpoint", 1, Inf, "#a445ee"),
  changeLabel(">0breakpoints", 1, Inf, "#a445ee"),
  changeLabel(">0changes", 1, Inf, "#a445ee"),
  changeLabel("normal", 0, 0, "#f6f4bf"),
  changeLabel("0breakpoints", 0, 0, "#f6f4bf"),
  changeLabel("0changes", 0, 0, "#f6f4bf"),
  changeLabel("1breakpoint", 1, 1, "#ff7d7d"),
  changeLabel("1change", 1, 1, "#ff7d7d"))
with(change.labels, stopifnot(min.changes <= max.changes))
stopifnot(0 <= change.labels$min.changes)
change.labels$possible.fn <- ifelse(0 < change.labels$min.changes, 1, 0)
change.labels$possible.fp <- ifelse(Inf == change.labels$max.changes, 0, 1)
setkey(change.labels, annotation)
change.colors <- paste(change.labels$color)
names(change.colors) <- change.labels$annotation

