\name{modelSelection}
\alias{modelSelection}
\title{Compute exact model selection function}
\description{Given loss.vec L_i, model.complexity K_i, the model selection
function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
the solutions (i, min.lambda, max.lambda) with i being the
solution for every lambda in (min.lambda, max.lambda). Use this
function after having computed changepoints and \code{loss} values for
each model, and before using \code{\link{labelError}}. This function uses the
linear time algorithm implemented in C code (\code{\link{modelSelectionC}}).}
\usage{modelSelection(models, loss = "loss", complexity = "complexity")}
\arguments{
  \item{models}{data.frame with one row per model. There must be at least two
columns models[[loss]] and models[[complexity]], but there can
also be other meta-data columns.}
  \item{loss}{character: column name of \code{models} to interpret as \code{loss} L_i.}
  \item{complexity}{character: column name of \code{models} to interpret as \code{complexity} K_i.}
}

\value{data.frame with a row for each model that can be selected for at
least one lambda value, and the following columns. (min.lambda,
max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
optimal penalty constants, on the original and log scale; the
other columns (and rownames) are taken from \code{models}. This should be
used as the \code{models} argument of \code{\link{labelError}}.}

\author{Toby Dylan Hocking}





