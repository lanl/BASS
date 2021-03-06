% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bassPCA.R
\name{sobolBasis}
\alias{sobolBasis}
\title{BASS Sensitivity Analysis}
\usage{
sobolBasis(
  mod,
  int.order,
  prior = NULL,
  mcmc.use = NULL,
  nind = NULL,
  n.cores = 1,
  parType = "fork",
  plot = F,
  verbose = T
)
}
\arguments{
\item{mod}{output from the \code{bassBasis} or \code{bassPCA} function.}

\item{int.order}{an integer indicating the highest order of interactions to include in the Sobol decomposition.}

\item{prior}{a list with the same number of elements as there are inputs to mod.  Each element specifies the prior for the particular input.  Each prior is specified as a list with elements \code{dist} (one of \code{c("normal", "student", "uniform")}), \code{trunc} (a vector of dimension 2 indicating the lower and upper truncation bounds, taken to be the data bounds if omitted), and for "normal" or "student" priors, \code{mean} (scalar mean of the Normal/Student, or a vector of means for a mixture of Normals or Students), \code{sd} (scalar standard deviation of the Normal/Student, or a vector of standard deviations for a mixture of Normals or Students), \code{df} (scalar degrees of freedom of the Student, or a vector of degrees of freedom for a mixture of Students), and \code{weights} (a vector of weights that sum to one for the mixture components, or the scalar 1).  If unspecified, a uniform is assumed with the same bounds as are represented in the input to mod.}

\item{mcmc.use}{an integer indicating which MCMC iteration to use for sensitivity analysis. Defaults to the last iteration.}

\item{nind}{number of Sobol indices to keep (will keep the largest nind).}

\item{n.cores}{number of cores to use (nearly linear speedup for adding cores).}

\item{parType}{either "fork" or "socket".  Forking is typically faster, but not compatible with Windows. If \code{n.cores==1}, \code{parType} is ignored.}

\item{plot}{logical; whether to plot results.}

\item{verbose}{logical; print progress.}
}
\value{
If non-functional (\code{func.var = NULL}), a list with two elements:
 \item{S}{a data frame of sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular main effect or interaction.  The values are the proportion of variance in the model that is due to each main effect or interaction.}
 \item{T}{a data frame of total sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular variable.}
 Otherwise, a list with four elements:
 \item{S}{an array with first dimension corresponding to MCMC samples (same length as \code{mcmc.use}), second dimension corresponding to different main effects and interactions (labeled in \code{names.ind}), and third dimension corresponding to the grid used for the functional variable.  The elements of the array are sensitivity indices.}
 \item{S.var}{same as \code{S}, but scaled in terms of total variance rather than percent of variance.}
 \item{names.ind}{a vector of names of the main effects and interactions used.}
}
\description{
Decomposes the variance of the BASS model into variance due to main effects, two way interactions, and so on, similar to the ANOVA decomposition for linear models.  Uses the Sobol' decomposition, which can be done analytically for MARS models.
}
\details{
Performs analytical Sobol' decomposition for each MCMC iteration in mcmc.use (each corresponds to a MARS model), yeilding a posterior distribution of sensitivity indices.  Can obtain Sobol' indices as a function of one functional variable.
}
\examples{
# See examples in bass documentation.
}
\seealso{
\link{bassPCA} and \link{bassBasis} for model fitting and \link{predict.bassBasis} for prediction.
}
\keyword{Sobol}
\keyword{decomposition}
