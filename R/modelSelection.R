changeLabel <- function(annotation, min.changes, max.changes, color){
  data.frame(annotation, min.changes, max.changes, color, row.names=annotation)
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
change.colors <- paste(change.labels$color)
names(change.colors) <- rownames(change.labels)

modelSelection <- structure(function # Exact model selection function
### Given a data.frame with "loss" column L_i, "model.complexity"
### column K_i, model selection function i*(lambda) = argmin_i L_i +
### lambda*K_i, compute all of the solutions (i, min.lambda,
### max.lambda) with i being the solution for every lambda in
### (min.lambda, max.lambda).
(){
  stopifnot(is.numeric(loss.vec))
  stopifnot(is.numeric(model.complexity))
  stopifnot(0 < diff(model.complexity))
  stopifnot(diff(loss.vec) < 0)
  stopifnot(length(loss.vec) == length(model.complexity))
  stopifnot(length(model.id) == length(model.complexity))
  n.models <- length(loss.vec)
  Kmax <- model.complexity[n.models]
  Kcurrent <- Kmax
  Lcurrent <- 0
  vK <- Kmax
  vL <- 0
  vP <- model.id[n.models]
  i <- 2
  min.complexity <- model.complexity[1]
  while(Kcurrent > min.complexity) {
    is.smaller <- model.complexity < Kcurrent
    is.current <- model.complexity == Kcurrent
    smallerK <- model.complexity[is.smaller]
    smallerID <- model.id[is.smaller]
    loss.term <- loss.vec[is.current] - loss.vec[is.smaller]
    complexity.term <- smallerK - model.complexity[is.current]
    lambdaTransition <- loss.term/complexity.term
    next.i <- which.min(lambdaTransition)
    Kcurrent <- smallerK[next.i]
    Lcurrent <- min(lambdaTransition)
    vL[i] <- Lcurrent
    vK[i] <- Kcurrent
    vP[i] <- smallerID[next.i]
    i <- i + 1
  }
  L <- log(vL)
  data.frame(
    min.lambda = vL,
    max.lambda = c(vL[-1], Inf),
    min.log.lambda = L,
    max.log.lambda = c(L[-1], Inf),
    model.complexity = vK,
    model.id=vP,
    model.loss=loss.vec[vP],
    row.names=vP)
### data.frame with a row for each model that can be selected for at
### least one lambda value, and the following columns. (min.lambda,
### max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
### optimal penalty constants, on the original and log scale;
### model.complexity are the K_i values; model.id are the model
### identifiers (also used for row names); and model.loss are the C_i
### values.
},ex=function(){

  if(require(neuroblastoma) && require(cghseg)){

    data(neuroblastoma, envir=environment())
    pro <- subset(neuroblastoma$profiles, profile.id==1 & chromosome=="X")
    max.segments <- 20
    fit <- cghseg:::segmeanCO(pro$logratio, Kmax=max.segments)
    seg.vec <- 1:max.segments
    exact.df <- modelSelection(fit$J.est, seg.vec, seg.vec)

    ## Solve the optimization using grid search.
    L.grid <- with(exact.df,{
      seq(min(max.log.lambda)-1,
          max(min.log.lambda)+1,
          l=100)
    })
    lambda.grid <- exp(L.grid)
    kstar.grid <- sapply(lambda.grid, function(lambda){
      crit <- with(exact.df, model.complexity * lambda + model.loss)
      picked <- which.min(crit)
      exact.df$model.id[picked]
    })
    grid.df <- data.frame(log.lambda=L.grid, segments=kstar.grid)

    if(require(ggplot2)){
      ## Compare the results.
      ggplot()+
        ggtitle("grid search (red) agrees with exact path computation (black)")+
        geom_segment(aes(min.log.lambda, model.id,
                         xend=max.log.lambda, yend=model.id),
                     data=exact.df)+
        geom_point(aes(log.lambda, segments),
                   data=grid.df, color="red", pch=1)+
        ylab("optimal model complexity (segments)")+
        xlab("log(lambda)")
    }
    
  }
  
})

largestContinuousMinimum <- structure(function
### Find the run of minimum cost with the largest size.
(cost,
 size
 ){
  m <- min(cost)
  is.min <- cost == m
  d <- c(diff(c(FALSE,is.min,FALSE)))
  ##print(data.frame(cost=c(cost,NA),size=c(size,NA),diff=d))
  starts <- which(d==1)
  ends <- which(d==-1)-1
  runs <- data.frame(starts,ends)
  ##print(runs)
  runs$size <- sapply(seq_along(starts),function(i){
    sum(size[ starts[i]:ends[i] ])
  })
  ##print(runs)
  largest <- which.max(runs$size)
  list(start=starts[largest],end=ends[largest])
}, ex=function(){

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

})

