#' Summary for fitted model ouput.
#'
#' @param object Fitted model object.
#' @param ... Other arguments.
#'
#'
#' @return Summary of fitted model object.
#' @export summary.fMmsm
#' @export
#'
#'
summary.fMmsm = function(object, ...){


  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL)

  Mmsm.post.object = object$Mmsm.post.object$Mmsm.post.object

  suStf1 = object$Mmsm.fit.object$fit$proc1$suStf
  suStf2 = object$Mmsm.fit.object$fit$proc2$suStf

  ii = 1
  jj = 1 # counter for smooths table

  # ************ #
  # PROCESS 1 ####
  # ************ #

  # COMPUTING P-VALUES OF COVARIATES THAT ENTER LINEARLY ####
  for(i in 1:sum(suStf1$whereQ != 0)){

    if(sum(suStf1$whereQ != 0) > 15) warning('The current code only displays results for up to 15 transitions. Contact maintainer for further details.')

    SE = sqrt(diag(Mmsm.post.object[[ii]]$Vp))

    this.trans.idx = suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)
    these.smooths = suStf1$start.pos.par.only.smooth >= suStf1$start.pos.par[i] & suStf1$start.pos.par.only.smooth < suStf1$start.pos.par[i+1] # identify smooths related to *this* transition - I think min is needed in case we have more than one smooth
    if(any(these.smooths)){
      not.smooths.idx = match(suStf1$start.pos.par[i]:(min(suStf1$start.pos.par.only.smooth[these.smooths])-1), this.trans.idx) # all linear terms will be before and all smooths terms after so we can select like this
    } else { # if there are no smooths then take all indices
      not.smooths.idx = match(this.trans.idx, this.trans.idx)
    }

    estimate = Mmsm.post.object[[ii]]$coefficients[not.smooths.idx]
    these.cov.names = Mmsm.post.object[[ii]]$coef.names[not.smooths.idx]
    se = SE[not.smooths.idx]
    ratio <- estimate/se
    pv <- 2 * pnorm(abs(ratio), lower.tail = FALSE)

    table[[ii]] = cbind(estimate, se, ratio, pv)
    dimnames(table[[ii]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    dimnames(table[[ii]])[[1]] = these.cov.names



    # COMPUTING P-VALUES OF SMOOTHS ####
    XX <- object$Mmsm.post.object$R

    edf.each.smooth = smooth.names = c()

    pTerms.df <- pTerms.chi.sq <- pTerms.pv <- list(0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0)


    this.trans.idx = suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)
    these.smooths = suStf1$start.pos.par.only.smooth >= suStf1$start.pos.par[i] & suStf1$start.pos.par.only.smooth < suStf1$start.pos.par[i+1] # identify smooths related to *this* transition

    if(any(these.smooths)){
      smooths.idx = match(min(suStf1$start.pos.par.only.smooth[these.smooths]):(suStf1$start.pos.par[i+1]-1), this.trans.idx) # all linear terms will be before and all smooths terms after so we can select like this - I think min is needed in case we have more than one smooth

      for(k in 1:sum(these.smooths)){

        this.smooth.start = suStf1$start.pos.par.only.smooth[these.smooths][k]
        this.smooth.end = ifelse(k == sum(these.smooths), max(this.trans.idx[smooths.idx]), suStf1$start.pos.par.only.smooth[these.smooths][k+1]-1)

        smooth.idx.each = match(this.smooth.start:this.smooth.end, this.trans.idx)

        edf.each.smooth[k] = sum(Mmsm.post.object[[ii]]$edf[smooth.idx.each])
        smooth.names[k] = Mmsm.post.object[[ii]]$smooth[[k]]$label

        if (max(Mmsm.post.object[[ii]]$smooth[[k]]$null.space.dim) == 0) { # is the max needed here? keep for now
          LRB <- rbind(XX, t(mroot(object$msm.fit.object$fit$S.h)))
          LRB <- cbind(LRB[, -smooth.idx.each], LRB[, smooth.idx.each])
          ind1 <- (ncol(LRB) - length(smooth.idx.each) + 1):ncol(LRB)
          Rm <- qr.R(qr(LRB, tol = 0, LAPACK = FALSE))[ind1,ind1]
          B <- mroot(Mmsm.post.object[[ii]]$Ve[smooth.idx.each, smooth.idx.each, drop = FALSE])
          b.hat <- Mmsm.post.object[[ii]]$coefficients[smooth.idx.each]
          d <- Rm %*% b.hat
          stat <- sum(d^2)
          ev <- eigen(crossprod(Rm %*% B), symmetric = TRUE,only.values = TRUE)$values
          ev[ev < 0] <- 0
          rank <- sum(ev > max(ev) * .Machine$double.eps^0.8)
          liu2 <- utils::getFromNamespace('liu2', 'mgcv')
          pval <- liu2(stat, ev)
          Tp <- list(stat = stat, pval = pval, rank = rank)
        }


        if (max(Mmsm.post.object[[ii]]$smooth[[k]]$null.space.dim) != 0) {
          b = Mmsm.post.object[[ii]]$coefficients[smooth.idx.each]
          V <- Mmsm.post.object[[ii]]$Vp[smooth.idx.each, smooth.idx.each, drop = FALSE]
          Xt <- XX[, smooth.idx.each, drop = FALSE]
          pTerms.df[[ii]][k] <- min(ncol(Xt), sum(Mmsm.post.object[[ii]]$edf[smooth.idx.each]))
          testStat <- utils::getFromNamespace('testStat', 'mgcv')
          Tp <- testStat(b, Xt, V, pTerms.df[[ii]][k], type = 0, res.df = -1)
        }

        pTerms.chi.sq[[ii]][k] <- Tp$stat
        pTerms.df[[ii]][k] <- Tp$rank
        pTerms.pv[[ii]][k] <- Tp$pval
      }

      tableN[[ii]] <- cbind(edf.each.smooth, pTerms.df[[ii]],
                            pTerms.chi.sq[[ii]], pTerms.pv[[ii]])
      dimnames(tableN[[ii]])[[2]] <- c("edf", "Ref.df",
                                       "Chi.sq", "p-value")
      dimnames(tableN[[ii]])[[1]] = smooth.names

    }
    # else {
    #   next # so the tableN entry just stays NULL as intialised above
    # }

    ii = ii + 1

  }





  # ************ #
  # PROCESS 2 ####
  # ************ #

  # COMPUTING P-VALUES OF COVARIATES THAT ENTER LINEARLY
  for(i in 1:sum(suStf2$whereQ != 0)){

    if(sum(suStf2$whereQ != 0) > 15) warning('The current code only displays results for up to 15 transitions. Contact maintainer for further details.')

    SE = sqrt(diag(Mmsm.post.object[[ii]]$Vp))

    this.trans.idx = suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)
    these.smooths = suStf2$start.pos.par.only.smooth >= suStf2$start.pos.par[i] & suStf2$start.pos.par.only.smooth < suStf2$start.pos.par[i+1] # identify smooths related to *this* transition - I think min is needed in case we have more than one smooth
    if(any(these.smooths)){
      not.smooths.idx = match(suStf2$start.pos.par[i]:(min(suStf2$start.pos.par.only.smooth[these.smooths])-1), this.trans.idx) # all linear terms will be before and all smooths terms after so we can select like this
    } else { # if there are no smooths then take all indices
      not.smooths.idx = match(this.trans.idx, this.trans.idx)
    }

    estimate = Mmsm.post.object[[ii]]$coefficients[not.smooths.idx]
    these.cov.names = Mmsm.post.object[[ii]]$coef.names[not.smooths.idx]
    se = SE[not.smooths.idx]
    ratio <- estimate/se
    pv <- 2 * pnorm(abs(ratio), lower.tail = FALSE)

    table[[ii]] = cbind(estimate, se, ratio, pv)
    dimnames(table[[ii]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    dimnames(table[[ii]])[[1]] = these.cov.names


    # COMPUTING P-VALUES OF SMOOTHS ####
    XX <- object$Mmsm.post.object$R


    edf.each.smooth = smooth.names = c()

    pTerms.df <- pTerms.chi.sq <- pTerms.pv <- list(0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0)


    this.trans.idx = suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)
    these.smooths = suStf2$start.pos.par.only.smooth >= suStf2$start.pos.par[i] & suStf2$start.pos.par.only.smooth < suStf2$start.pos.par[i+1] # identify smooths related to *this* transition

    if(any(these.smooths)){
      smooths.idx = match(min(suStf2$start.pos.par.only.smooth[these.smooths]):(suStf2$start.pos.par[i+1]-1), this.trans.idx) # all linear terms will be before and all smooths terms after so we can select like this - I think min is needed in case we have more than one smooth

      for(k in 1:sum(these.smooths)){

        this.smooth.start = suStf2$start.pos.par.only.smooth[these.smooths][k]
        this.smooth.end = ifelse(k == sum(these.smooths), max(this.trans.idx[smooths.idx]), suStf2$start.pos.par.only.smooth[these.smooths][k+1]-1)

        smooth.idx.each = match(this.smooth.start:this.smooth.end, this.trans.idx)

        edf.each.smooth[k] = sum(Mmsm.post.object[[ii]]$edf[smooth.idx.each])
        smooth.names[k] = Mmsm.post.object[[ii]]$smooth[[k]]$label

        if (max(Mmsm.post.object[[ii]]$smooth[[k]]$null.space.dim) == 0) { # is the max needed here? keep for now
          LRB <- rbind(XX, t(mroot(object$msm.fit.object$fit$S.h)))
          LRB <- cbind(LRB[, -smooth.idx.each], LRB[, smooth.idx.each])
          ind1 <- (ncol(LRB) - length(smooth.idx.each) + 1):ncol(LRB)
          Rm <- qr.R(qr(LRB, tol = 0, LAPACK = FALSE))[ind1,ind1]
          B <- mroot(Mmsm.post.object[[ii]]$Ve[smooth.idx.each, smooth.idx.each, drop = FALSE])
          b.hat <- Mmsm.post.object[[ii]]$coefficients[smooth.idx.each]
          d <- Rm %*% b.hat
          stat <- sum(d^2)
          ev <- eigen(crossprod(Rm %*% B), symmetric = TRUE,only.values = TRUE)$values
          ev[ev < 0] <- 0
          rank <- sum(ev > max(ev) * .Machine$double.eps^0.8)
          liu2 <- utils::getFromNamespace('liu2', 'mgcv')
          pval <- liu2(stat, ev)
          Tp <- list(stat = stat, pval = pval, rank = rank)
        }


        if (max(Mmsm.post.object[[ii]]$smooth[[k]]$null.space.dim) != 0) {
          b = Mmsm.post.object[[ii]]$coefficients[smooth.idx.each]
          V <- Mmsm.post.object[[ii]]$Vp[smooth.idx.each, smooth.idx.each, drop = FALSE]
          Xt <- XX[, smooth.idx.each, drop = FALSE]
          pTerms.df[[ii]][k] <- min(ncol(Xt), sum(Mmsm.post.object[[ii]]$edf[smooth.idx.each]))
          testStat <- utils::getFromNamespace('testStat', 'mgcv')
          Tp <- testStat(b, Xt, V, pTerms.df[[ii]][k], type = 0, res.df = -1)
        }

        pTerms.chi.sq[[ii]][k] <- Tp$stat
        pTerms.df[[ii]][k] <- Tp$rank
        pTerms.pv[[ii]][k] <- Tp$pval
      }

      tableN[[ii]] <- cbind(edf.each.smooth, pTerms.df[[ii]],
                            pTerms.chi.sq[[ii]], pTerms.pv[[ii]])
      dimnames(tableN[[ii]])[[2]] <- c("edf", "Ref.df",
                                       "Chi.sq", "p-value")
      dimnames(tableN[[ii]])[[1]] = smooth.names

    }
    # else {
    #   next # so the tableN entry just stays NULL as intialised above
    # }


    ii = ii + 1
  }



  state.pairs.CT.out1 = state.pairs.CT(data = suStf1$data, whereQ = suStf1$whereQ, nstates = suStf1$nstates,
                                      time = suStf1$tte, state = '(state)', id = '(id)')
  state.pairs.CT.out2 = state.pairs.CT(data = suStf2$data, whereQ = suStf2$whereQ, nstates = suStf2$nstates,
                                      time = suStf2$tte, state = '(state)', id = '(id)')


  idx.lab1 = which(t(suStf1$whereQ) != 0, arr.ind = T)
  idx.lab2 = which(t(suStf2$whereQ) != 0, arr.ind = T)


  # REVIEW OUTPUT TO ONLY KEEP THE THINGS WE ACTUALLY HAVE
  res <- list(table = table, tableN = tableN,
              tableP1 = table[[1]], tableP2 = table[[2]], tableP3 = table[[3]], tableP4 = table[[4]], tableP5 = table[[5]],
              tableP6 = table[[6]], tableP7 = table[[7]], tableP8 = table[[8]], tableP9 = table[[9]], tableP10 = table[[10]],
              tableP11 = table[[11]], tableP12 = table[[12]], tableP13 = table[[13]], tableP14 = table[[14]], tableP15 = table[[15]],
              tableNP1 = tableN[[1]], tableNP2 = tableN[[2]], tableNP3 = tableN[[3]], tableNP4 = tableN[[4]], tableNP5 = tableN[[5]],
              tableNP6 = tableN[[6]], tableNP7 = tableN[[7]], tableNP8 = tableN[[8]], tableNP9 = tableN[[9]], tableNP10 = tableN[[10]],
              tableNP11 = tableN[[11]], tableNP12 = tableN[[12]], tableNP13 = tableN[[13]], tableNP14 = tableN[[14]], tableNP15 = tableN[[15]],
              n1 = object$Mmsm.fit.object$fit$proc1$n, n2 = object$Mmsm.fit.object$fit$proc2$n,
              n = object$Mmsm.fit.object$fit$proc1$n + object$Mmsm.fit.object$fit$proc2$n,
              N = object$Mmsm.fit.object$fit$proc1$N,
              # sigma2 = object$sigma2, sigma = object$sigma2, Model = object$Model,
              # nu = object$nu, sigma2.a = object$sigma2.a, sigma.a = object$sigma2.a,
              # nu.a = object$nu.a,
              idx.lab1 = idx.lab1, idx.lab2 = idx.lab2,
              short.formula1 = object$Mmsm.fit.object$fit$proc1$short.formula, short.formula2 = object$Mmsm.fit.object$fit$proc2$short.formula,
              suStf1 = suStf1, suStf2 = suStf2,
              # formula1 = formula1, formula2 = formula2, formula3 = formula3, formula4 = formula4, formula5 = formula5,
              # formula6 = formula6, formula7 = formula7, formula8 = formula8, formula9 = formula9, formula10 = formula10,
              # formula11 = formula11, formula12 = formula12, formula13 = formula13, formula14 = formula14, formula15 = formula15,
              t.edf = object$Mmsm.post.object$t.edf,
              t.edf1 = object$Mmsm.post.object$t.edf1, t.edf2 = object$Mmsm.post.object$t.edf2,
              postpost = object$Mmsm.post.object$Mmsm.post.object,
              state.pairs.CT.out1 = state.pairs.CT.out1, state.pairs.CT.out2 = state.pairs.CT.out2,
              phiHat = exp(object$Mmsm.fit.object$fit$argument[length(object$Mmsm.fit.object$fit$argument)])
              # CImu = CImu, mu = mu, CIsig = CIsig2, CInu = CInu, margins = object$margins,
              # l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3,
              # l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6,
              # l.sp7 = object$l.sp7, l.sp8 = object$l.sp8, X2.null = is.null(object$X2),
              # univar.gamlss = TRUE, surv.flex = object$surv.flex,
              # K1 = NULL, robust = object$robust, indx = object$fit$indx
  )

  class(res) <- "summary.fMmsm"
  res

}
