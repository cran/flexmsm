#' Convergence diagnostics on fitted model output.
#'
#' @param object Fitted model object.
#' @param also.unpen If \code{TRUE}, displays eigenvalues range also for the unpenalised Hessian. Defaults to \code{FALSE}.
#'
#'
#' @return Convergence diagnostics.
#' @export
#'
#'


conv.check = function(object, also.unpen = FALSE){

  fit = object$msm.fit.object$fit

  e.v <- eigen(fit$hessian, symmetric = TRUE, only.values = TRUE)$values
  e.v.np <- eigen(fit$hessian - fit$S.h, symmetric = TRUE, only.values = TRUE)$values

  cat("\nLargest absolute gradient value:", max(abs(fit$gradient)), '\n')


  if (min(e.v) > 0) {
    cat("\nPenalized analytic Hessian is positive definite\n",
        sep = "")
  } else {
    cat("\nPenalized analytic Hessian is not positive definite\n",
        sep = "")
  }
  cat("Eigenvalue range: [", min(e.v), ",", max(e.v), "]\n", sep = "")

  if(also.unpen){ # display eigenvalues information also for unpenalised hessian
    if (min(e.v.np) > 0) {
      cat("\nUnpenalized analytic Hessian is positive definite\n",
          sep = "")
    } else {
      cat("\nUnpenalized analytic Hessian is not positive definite\n",
          sep = "")
    }
    cat("Eigenvalue range (unpenalized): [", min(e.v.np), ",", max(e.v.np), "]\n", sep = "")
  }


  if(length(object$msm.fit.object$Qmatr.diagnostics.list) > 0) cat("\nQ matrix range: [", min(object$msm.fit.object$Qmatr.diagnostics.list[[length(object$msm.fit.object$Qmatr.diagnostics.list)]]), ",",
                                                                   max(object$msm.fit.object$Qmatr.diagnostics.list[[length(object$msm.fit.object$Qmatr.diagnostics.list)]]), "]\n", sep = "")

  cat("\nTrust region iterations before smoothing parameter estimation:", object$msm.fit.object$fit.all[[1]]$iterations)
  cat("\nLoops for smoothing parameter estimation:", object$msm.fit.object$iter.sp)
  cat("\nTrust region iterations within smoothing loops:", object$msm.fit.object$iter.inner)


  cat("\n\n")

}
