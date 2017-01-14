\name{modelSelection}
\alias{modelSelection}
\title{modelSelection}
\description{Given loss.vec L_i, model.complexity K_i, the model selection
function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
the solutions (i, min.lambda, max.lambda) with i being the
solution for every lambda in (min.lambda, max.lambda). This
function uses the quadratic time algorithm implemented in R code.}
\usage{modelSelection(models, loss = "loss", complexity = "complexity")}
\arguments{
  \item{models}{data.frame with one row per model. There must be at
least two columns [[loss]] and [[complexity]], but there can
also be other meta-data columns.}
  \item{loss}{character: column name of models to interpret as loss L_i.}
  \item{complexity}{character: column name of models to interpret as complexity K_i.}
}

\value{data.frame with a row for each model that can be selected for at
least one lambda value, and the following columns. (min.lambda,
max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
optimal penalty constants, on the original and log scale;
the other columns (and rownames) are taken from models.}

\author{Toby Dylan Hocking}





