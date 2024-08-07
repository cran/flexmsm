% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmsmObject.R
\name{fmsmObject}
\alias{fmsmObject}
\title{Fitted fmsmObject object}
\value{
\item{suStf}{A list with all of the quantities used for estimation and post-estimation computations. This includes the full design matrix \code{full.X}, the starting parameters used, \code{params.0} and \code{sp.0}, and more technical quantities such as the positions of the smooths' parameters and of the parametric coefficients.}
\item{msm.fit.object}{This contains all of the details of the model fitting.}
\item{msm.post.object}{This contains all of the post-estimation details.}
\item{formula}{Formula used in the model specification.}
\item{short.formula}{Short version of the model specification, i.e. only non-zero transition specifications are included.}
\item{n}{Number of observations in the dataset.}
\item{N}{Number of unique individuals.}
\item{logLik}{The value of the log-likelihood at convergence.}
\item{t.edf}{Total effective degrees of freedom.}
\item{singleComp}{If present, log-likelihood, gradient and Hessian computed at the starting parameter (without fitting).}
}
\description{
The \code{\link{fmsm}} function returns the fitted model object \code{fmsmObject}. This is of class "fmsm" and includes the components listed below. These are intended for confident users. To extract results
from the fitted model objects, functions such as \code{\link{summary.fmsm}}, \code{\link{plot.fmsm}}, \code{\link{Q.pred}} and \code{\link{P.pred}} should be used instead.
}
\seealso{
\code{\link{fmsm}}, \code{\link{summary.fmsm}}
}
