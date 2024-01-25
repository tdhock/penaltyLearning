labelError <- structure(function # Compute incorrect labels
### Compute incorrect labels for several change-point detection
### problems and models. Use this function after having computed
### changepoints, loss values, and model selection functions
### (see modelSelection). The next step after labelError is typically
### computing target intervals of log(penalty) values that predict
### changepoints with minimum incorrect labels for each problem (see
### targetIntervals).
(models,
### data.frame with one row per (problem,model) combination, typically
### the output of modelSelection(...). There is a row for each
### changepoint model that could be selected for a particular
### segmentation problem. There should be columns problem.vars (for
### problem ID) and model.vars (for model complexity).
  labels,
### data.frame with one row per (problem,region). Each label defines a
### region in a particular segmentation problem, and a range of
### predicted changepoints which are consistent in that region. There
### should be a column "annotation" with takes one of the
### corresponding values in the annotation column of change.labels
### (used to determine the range of predicted changepoints which are
### consistent). There should also be a columns problem.vars (for
### problem ID) and label.vars (for region start/end).
  changes,
### data.frame with one row per (problem,model,change), for each
### predicted changepoint (in each model and segmentation
### problem). Should have columns problem.vars (for problem ID),
### model.vars (for model complexity), and change.var (for changepoint
### position).
  change.var="chromStart",
### character(length=1): column name of predicted change-point
### position in labels. The default "chromStart" is useful for genomic
### data with segment start/end positions stored in columns named
### chromStart/chromEnd. A predicted changepoint at position X is
### interpreted to mean a changepoint between X and X+1.
  label.vars=c("min", "max"),
### character(length=2): column names of start and end positions of
### labels, in same units as change-point positions. The default is
### c("min", "max"). Labeled regions are (start,end] -- open on the
### left and closed on the right, so for example a 0changes annotation
### between start=10 and end=20 means that any predicted changepoint
### at 11, ..., 20 is a false positive.
  model.vars="n.segments",
### character: column names used to identify model complexity. The
### default "n.segments" is for change-point models such as in the
### jointseg and changepoint packages.
  problem.vars=character(0),
### character: column names used to identify data set / segmentation
### problem, should be present in all three data tables (models,
### labels, changes).
  annotations=change.labels
### data.table with columns annotation, min.changes, max.changes,
### possible.fn, possible.fp which is joined to labels in order to
### determine how to compute false positives and false negatives for
### each annotation.
) {
  weight <- annotation <- fp <- max.changes <- pred.changes <- fn <-
    min.changes <- status <- possible.fp <- possible.fn <- NULL
### The code above is to avoid CRAN NOTEs like
### labelError: no visible binding for global variable
  stopifnot(is.data.frame(models))
  stopifnot(is.data.frame(labels))
  stopifnot(is.data.frame(changes))
  if(!(
    is.character(problem.vars) &&
      0 < length(problem.vars) &&
      all(!is.na(problem.vars)) &&
      all(problem.vars %in% names(changes)) &&
      all(problem.vars %in% names(labels)) &&
      all(problem.vars %in% names(models))
  )){
    stop("problem.vars should be a character vector of column names present in models, changes, and labels (ID for separate changepoint detection problems)")
  }
  if(!(
    is.character(label.vars) &&
    length(label.vars)==2 &&
    all(label.vars %in% names(labels))
  )){
    stop("label.vars should be a 2-element character vector of labels column names (start and end of labeled region)")
  }
  if(any(labels[[ label.vars[[2]] ]] <= labels[[ label.vars[[1]] ]])){
    stop("label start must be less than end")
  }
  if(!(
    is.character(change.var) &&
    length(change.var)==1 &&
    change.var %in% names(changes)
  )){
    stop("change.var should be a column name of changes (position of predicted changepoints)")
  }
  if(!(
    is.character(model.vars) &&
      0 < length(model.vars) &&
      all(model.vars %in% names(models)) &&
      all(model.vars %in% names(changes))
  )){
    stop("model.vars should be a column name of both models and changes (ID for model complexity, typically the number of changepoints or segments)")
  }
  labels.dt <- data.table(labels)
  setkeyv(labels.dt, c(problem.vars, label.vars))
  labels.dt[, {
    end <- .SD[[ label.vars[[2]] ]][-.N]
    next.start <- .SD[[ label.vars[[1]] ]][-1]
    if(any(next.start < end)){
      stop("each label end must be <= next label start")
    }
  }, by=problem.vars]
  if("weight" %in% names(labels.dt)){
    stopifnot(is.numeric(labels.dt$weight))
    stopifnot(0 < labels.dt$weight)
  }else{
    labels.dt[, weight := 1 ]
  }
  labels.info <- annotations[labels.dt, on=list(annotation), nomatch=0L]
  if(nrow(labels.info)!=nrow(labels.dt)){
    stop("labels$annotation must be one of annotations$annotation")
  }
  models.dt <- data.table(models)
  model.labels <- models.dt[labels.info, on=problem.vars, allow.cartesian=TRUE]
  if(any(is.na(model.labels[, model.vars, with=FALSE]))){
    stop("some labels have no models")
  }
  changes.dt <- data.table(changes)
  changes.per.problem <- changes.dt[models.dt, list(
    pred.changes=.N
  ), by=.EACHI, on=c(problem.vars, model.vars)]
  model.counts <- changes.per.problem[, list(
    models=.N
  ), by=c(problem.vars, "pred.changes")]
  bad.models <- model.counts[1 < models]
  if(nrow(bad.models)){
    print(bad.models)
    stop("each model should have a different number of changes, problems displayed above")
  }
  if(nrow(changes.dt)==0){
    over.dt <- data.table(
      model.labels[, c(problem.vars, model.vars, label.vars), with=FALSE],
      pred.changes=0L)
  }else{
    over.dt <- changes.dt[model.labels, list(
      pred.changes=.N
      ), by=.EACHI, on=c(
      problem.vars, model.vars,
      paste0(change.var, c(">", "<="), label.vars))]
    ## Is this a bug in data.table? Why should I have to set names back
    ## to start and end (they are both pos after the
    ## join). https://github.com/Rdatatable/data.table/issues/1700
    setnames(over.dt, c(
      problem.vars, model.vars,
      label.vars,
      "pred.changes"))
  }
  changes.per.label <- over.dt[model.labels, on=c(
    problem.vars, model.vars, label.vars)]
  changes.per.label[, fp := ifelse(max.changes < pred.changes, weight, 0)]
  changes.per.label[, fn := ifelse(pred.changes < min.changes, weight, 0)]
  changes.per.label[, status := ifelse(
    fp, "false positive", ifelse(
      fn, "false negative", "correct"))]
  error.totals <- changes.per.label[, list(
    possible.fp=sum(possible.fp*weight),
    fp=sum(fp),
    possible.fn=sum(possible.fn*weight),
    fn=sum(fn),
    labels=sum(weight),
    errors=sum(fp+fn)),
    by=c(problem.vars, model.vars)]
  list(
    model.errors=models.dt[error.totals, on=c(problem.vars, model.vars)],
    label.errors=changes.per.label)
### list of two data.tables: label.errors has one row for every
### combination of models and labels, with status column that
### indicates whether or not that model commits an error in that
### particular label; model.errors has one row per model, with columns
### for computing target intervals and ROC curves (see targetIntervals
### and ROChange).
}, ex=function() {

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

})

### Describe an annotated region label for supervised change-point detection.
changeLabel <- function(annotation, min.changes, max.changes, color){
  data.table(annotation, min.changes, max.changes, color)
}

### data.table of meta-data for label types.
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

### character vector of change-point label colors, to be used with
### ggplot2::scale_*_manual
change.colors <- paste(change.labels$color)
names(change.colors) <- change.labels$annotation

