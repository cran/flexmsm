#' Fitted fmsmObject object
#'
#' @name fmsmObject
#'
#' @description
#' The \code{\link{fmsm}} function returns the fitted model object \code{fmsmObject}. This is of class "fmsm" and includes the components listed below. These are intended for confident users. To extract results
#' from the fitted model objects, functions such as \code{\link{summary.fmsm}}, \code{\link{plot.fmsm}}, \code{\link{Q.pred}} and \code{\link{P.pred}} should be used instead.
#'
#' @return \tabular{ll}{
#' \code{suStf} \tab A list with all of the quantities used for estimation and post-estimation computations. This includes the full design matrix \code{full.X}, the starting parameters used, \code{params.0} and \code{sp.0}, and more technical quantities such as the positions of the smooths' parameters and of the parametric coefficients. \cr
#' \tab \cr
#' \code{msm.fit.object} \tab This contains all of the details of the model fitting. \cr
#' \tab \cr
#' \code{msm.post.object} \tab This contains all of the post-estimation details. \cr
#' \tab \cr
#' \code{formula} \tab Formula used in the model specification. \cr
#' \tab \cr
#' \code{short.formula} \tab Short version of the model specification, i.e. only non-zero transition specifications are included. \cr
#' \tab \cr
#' \code{n} \tab Number of observations in the dataset. \cr
#' \tab \cr
#' \code{N} \tab Number of unique individuals. \cr
#' \tab \cr
#' \code{logLik} \tab The value of the log-likelihood at convergence. \cr
#' \tab \cr
#' \code{t.edf} \tab Total effective degrees of freedom. \cr
#' }
#'
#'
#' @seealso \code{\link{fmsm}}, \code{\link{summary.fmsm}}
#'
#'
#'
NULL