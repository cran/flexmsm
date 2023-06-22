


#' Extract the log likelihood for the fitted multistate model
#'
#' @description
#' It extracts the log-likelihood for a fitted \code{fmsm} model.
#'
#'
#' @param object Fitted model object of class \code{fmsm} produced by function \code{\link{fmsm}}.
#' @param ... Unused in this case.
#'
#' @return Standard logLik object.
#'
#' @export logLik.fmsm
#' @export
#'


logLik.fmsm = function(object, ...){

  if (length(list(...))) {
    warning("extra arguments discarded")
  }

  lk <- object$logLik

  attr(lk, "nobs") <- object$n
  attr(lk, "df") <- object$t.edf
  class(lk) <- c("logLik")
  lk

}
