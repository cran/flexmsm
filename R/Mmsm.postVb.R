




Mmsm.postVb = function(Mmsm.fit.object){


  tolH <- sqrt(.Machine$double.eps)
  He <- HeSh <- Mmsm.fit.object$fit$hessian
  He.eig <- eigen(He, symmetric = TRUE)
  if (min(He.eig$values) < tolH && sign(min(sign(He.eig$values))) == -1) He.eig$values <- abs(He.eig$values)

  if (min(He.eig$values) < tolH) {
    pep <- which(He.eig$values < tolH)
    He.eig$values[pep] <- tolH
  }

  Vb <- He.eig$vectors %*% tcrossprod(diag(1/He.eig$values,
                                           nrow = length(He.eig$values), ncol = length(He.eig$values)),
                                      He.eig$vectors)

  Vb <- Vb1 <- (Vb + t(Vb))/2

  if (length(Mmsm.fit.object$sp) > 0) {

    HeSh <- He - Mmsm.fit.object$fit$S.h

    if (Mmsm.fit.object$sp.method == "perf") {
      F <- diag(Mmsm.fit.object$magpp$edf)
      F1 <- diag(Mmsm.fit.object$magpp$edf1)
      R <- Mmsm.fit.object$bs.mgfit$R
      Ve <- Vb %*% HeSh %*% Vb
      Ve <- (Ve + t(Ve))/2
    }


    if (Mmsm.fit.object$sp.method == "efs") {
      lbb <- Sl.initial.repara(Mmsm.fit.object$Sl, HeSh, inverse = FALSE)
      p <- ncol(lbb)
      ipiv <- piv <- attr(Mmsm.fit.object$L, "pivot")
      ipiv[piv] <- 1:p
      lbb <- Mmsm.fit.object$D * t(Mmsm.fit.object$D * lbb)
      R <- suppressWarnings(chol(lbb, pivot = TRUE))
      if (attr(R, "rank") < ncol(R)) {
        retry <- TRUE
        tol <- 0
        eh <- eigen(lbb, symmetric = TRUE)
        mev <- max(eh$values)
        dtol <- 1e-07
        while (retry) {
          eh$values[eh$values < tol * mev] <- tol *
            mev
          R <- sqrt(eh$values) * t(eh$vectors)
          lbb <- crossprod(R)
          Hp <- lbb + Mmsm.fit.object$D * t(Mmsm.fit.object$D *
                                             Mmsm.fit.object$St)
          Mmsm.fit.object$L <- suppressWarnings(chol(Hp,
                                                    pivot = TRUE))
          if (attr(Mmsm.fit.object$L, "rank") == ncol(Hp)) {
            R <- t(t(R)/Mmsm.fit.object$D)
            retry <- FALSE
          }
          else {
            tol <- tol + dtol
            dtol <- dtol * 10
          }
        }
      }
      else {
        ipiv <- piv <- attr(R, "pivot")
        ipiv[piv] <- 1:p
        R <- t(t(R[, ipiv])/Mmsm.fit.object$D)
      }
      R <- Sl.repara(Mmsm.fit.object$rp, R, inverse = TRUE,
                     both.sides = FALSE)
      R <- Sl.initial.repara(Mmsm.fit.object$Sl, R, inverse = TRUE,
                             both.sides = FALSE, cov = FALSE)
      F <- Vb %*% HeSh
      F1 <- diag(2 * diag(F) - rowSums(t(F) * F))
      Ve <- F %*% Vb
      Ve <- (Ve + t(Ve))/2
    }


  } else {
    Ve <- Vb
    F <- F1 <- diag(rep(1, dim(Vb)[1]))
    R <- Mmsm.fit.object$bs.mgfit$R
  }




  t.edf <- sum(diag(F))
  Vb11 <- F11 <- Mmsm.fit.object$fit$hessian
  d1 <- dim(Vb11)[1]
  d2 <- dim(Vb11)[2]
  Vb11[1:d1, 1:d2] <- Vb[1:d1, 1:d2]
  F11[1:d1, 1:d2] <- F[1:d1, 1:d2]
  Vb <- Vb11
  F <- F11


  dimnames(Mmsm.fit.object$fit$hessian)[[1]] <- dimnames(Mmsm.fit.object$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(Mmsm.fit.object$fit$argument)



  list(Vb1 = Vb1, He = He, Vb = Vb, HeSh = HeSh,
       F = F, F1 = F1, R = R, Ve = Ve, t.edf = t.edf,
       Mmsm.fit.object = Mmsm.fit.object)


}
