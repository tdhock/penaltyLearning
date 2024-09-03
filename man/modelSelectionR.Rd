\name{modelSelectionR}
\alias{modelSelectionR}
\title{Exact model selection function}
\description{Given \code{loss.vec} L_i, \code{model.complexity} K_i, the model selection
function i*(lambda) = argmin_i L_i + lambda*K_i, compute all of
the solutions (i, min.lambda, max.lambda) with i being the
solution for every lambda in (min.lambda, max.lambda). This
function uses the quadratic time algorithm implemented in R code.
This function is mostly meant for internal use and comparison --
it is instead recommended to use \code{\link{modelSelection}}.}
\usage{modelSelectionR(loss.vec, 
    model.complexity, 
    model.id)}
\arguments{
  \item{loss.vec}{numeric vector: loss L_i}
  \item{model.complexity}{numeric vector: model complexity K_i}
  \item{model.id}{vector: indices i}
}

\value{data.frame with a row for each model that can be selected for at
least one lambda value, and the following columns. (min.lambda,
max.lambda) and (min.log.lambda, max.log.lambda) are intervals of
optimal penalty constants, on the original and log scale;
\code{model.complexity} are the K_i values; \code{model.id} are the model
identifiers (also used for row names); and model.loss are the C_i
values.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

loss.vec <- c(
  -9.9, -12.8, -19.2, -22.1, -24.5, -26.1, -28.5, -30.1, -32.2, 
  -33.7, -35.2, -36.8, -38.2, -39.5, -40.7, -41.8, -42.8, -43.9, 
  -44.9, -45.8)
seg.vec <- seq_along(loss.vec)
penaltyLearning::modelSelectionR(loss.vec, seg.vec, seg.vec)

}
