#' Summary for fitted model ouput.
#'
#' @param object Fitted model object.
#' @param ... Other arguments.
#'
#'
#' @return Summary of fitted model object.
#' @export summary.fmsm
#' @export
#'
#'
summary.fmsm = function(object, ...){


  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL) # for now we limit ourselves to 15 transitions (never reached anyway)

  # formula1 = formula2 = formula3 = formula4 = formula5 = formula6 = formula7 = formula8 = formula9 = formula10 = formula11 = formula12 = formula13 = formula14 = formula15 = NULL

  idx.lab = which(t(object$suStf$whereQ) != 0, arr.ind = T)

  # # FOR NOW WE WILL IGNORE THE COMPUTATION OF THESE CONFIDENCE INTERVALS HERE INSIDE THE SUMMARY FUNCTION - perhaps change later
  # bs <- rMVN(n.sim, mean = object$coefficients, sigma = Vb)
  # susutsnR <- susutsn(object, bs, lf, cont1par, cont2par,
  #                     cont3par, prob.lev, type = "gamls")
  # CIsig2 <- susutsnR$CIsig21
  # CInu <- susutsnR$CInu1
  # CImu <- susutsnR$CImu
  # mu <- susutsnR$mu ***************************************************************************************************************

  msm.post.object = object$msm.post.object$msm.post.object
  suStf = object$suStf

  # COMPUTING P-VALUES OF COVARIATES THAT ENTER LINEARLY ####
  for(i in 1:sum(suStf$whereQ != 0)){

    if(sum(suStf$whereQ != 0) > 15) warning('The current code only displays results for up to 15 transitions. Contact maintainer for further details.')

    SE = sqrt(diag(msm.post.object[[i]]$Vp))

    this.trans.idx = suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)
    these.smooths = suStf$start.pos.par.only.smooth >= suStf$start.pos.par[i] & suStf$start.pos.par.only.smooth < suStf$start.pos.par[i+1] # identify smooths related to *this* transition - I think min is needed in case we have more than one smooth
    if(any(these.smooths)){
      not.smooths.idx = match(suStf$start.pos.par[i]:(min(suStf$start.pos.par.only.smooth[these.smooths])-1), this.trans.idx) # all linear terms will be before and all smooths terms after so we can select like this
    } else { # if there are no smooths then take all indices
      not.smooths.idx = match(this.trans.idx, this.trans.idx)
    }

    estimate = msm.post.object[[i]]$coefficients[not.smooths.idx]
    these.cov.names = msm.post.object[[i]]$coef.names[not.smooths.idx]
    se = SE[not.smooths.idx]
    ratio <- estimate/se
    pv <- 2 * pnorm(abs(ratio), lower.tail = FALSE)

    table[[i]] = cbind(estimate, se, ratio, pv)
    dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    dimnames(table[[i]])[[1]] = these.cov.names

  }



  # COMPUTING P-VALUES OF SMOOTHS ####
  XX <- object$msm.post.object$R

  for(i in 1:sum(suStf$whereQ != 0)){

    edf.each.smooth = smooth.names = c()

    pTerms.df <- pTerms.chi.sq <- pTerms.pv <- list(0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0)


    this.trans.idx = suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)
    these.smooths = suStf$start.pos.par.only.smooth >= suStf$start.pos.par[i] & suStf$start.pos.par.only.smooth < suStf$start.pos.par[i+1] # identify smooths related to *this* transition

    if(any(these.smooths)){
      smooths.idx = match(min(suStf$start.pos.par.only.smooth[these.smooths]):(suStf$start.pos.par[i+1]-1), this.trans.idx) # all linear terms will be before and all smooths terms after so we can select like this - I think min is needed in case we have more than one smooth

      for(k in 1:sum(these.smooths)){

        this.smooth.start = suStf$start.pos.par.only.smooth[these.smooths][k]
        this.smooth.end = ifelse(k == sum(these.smooths), max(this.trans.idx[smooths.idx]), suStf$start.pos.par.only.smooth[these.smooths][k+1]-1)

        smooth.idx.each = match(this.smooth.start:this.smooth.end, this.trans.idx)

        edf.each.smooth[k] = sum(msm.post.object[[i]]$edf[smooth.idx.each])
        smooth.names[k] = msm.post.object[[i]]$smooth[[k]]$label

        if (max(msm.post.object[[i]]$smooth[[k]]$null.space.dim) == 0) { # is the max needed here? keep for now
          LRB <- rbind(XX, t(mroot(object$msm.fit.object$fit$S.h)))
          LRB <- cbind(LRB[, -smooth.idx.each], LRB[, smooth.idx.each])
          ind1 <- (ncol(LRB) - length(smooth.idx.each) + 1):ncol(LRB)
          Rm <- qr.R(qr(LRB, tol = 0, LAPACK = FALSE))[ind1,ind1]
          B <- mroot(msm.post.object[[i]]$Ve[smooth.idx.each, smooth.idx.each, drop = FALSE])
          b.hat <- msm.post.object[[i]]$coefficients[smooth.idx.each]
          d <- Rm %*% b.hat
          stat <- sum(d^2)
          ev <- eigen(crossprod(Rm %*% B), symmetric = TRUE,only.values = TRUE)$values
          ev[ev < 0] <- 0
          rank <- sum(ev > max(ev) * .Machine$double.eps^0.8)
          liu2 <- utils::getFromNamespace('liu2', 'mgcv')
          pval <- liu2(stat, ev)
          Tp <- list(stat = stat, pval = pval, rank = rank)
        }


        if (max(msm.post.object[[i]]$smooth[[k]]$null.space.dim) != 0) {
          b = msm.post.object[[i]]$coefficients[smooth.idx.each]
          V <- msm.post.object[[i]]$Vp[smooth.idx.each, smooth.idx.each, drop = FALSE]
          Xt <- XX[, smooth.idx.each, drop = FALSE]
          pTerms.df[[i]][k] <- min(ncol(Xt), sum(msm.post.object[[i]]$edf[smooth.idx.each]))
          testStat <- utils::getFromNamespace('testStat', 'mgcv')
          Tp <- testStat(b, Xt, V, pTerms.df[[i]][k], type = 0, res.df = -1)
        }

        pTerms.chi.sq[[i]][k] <- Tp$stat
        pTerms.df[[i]][k] <- Tp$rank
        pTerms.pv[[i]][k] <- Tp$pval
      }

    } else {
      next # so the tableN entry just stays NULL as intialised above
    }





    tableN[[i]] <- cbind(edf.each.smooth, pTerms.df[[i]],
                    pTerms.chi.sq[[i]], pTerms.pv[[i]])
    dimnames(tableN[[i]])[[2]] <- c("edf", "Ref.df",
                               "Chi.sq", "p-value")
    dimnames(tableN[[i]])[[1]] = smooth.names
  }


  state.pairs.CT.out = state.pairs.CT(data = suStf$data, whereQ = suStf$whereQ, nstates = suStf$nstates,
                                      time = suStf$tte, state = '(state)', id = '(id)')



    # REVIEW OUTPUT TO ONLY KEEP THE THINGS WE ACTUALLY HAVE
    res <- list(table = table, tableN = tableN,
                tableP1 = table[[1]], tableP2 = table[[2]], tableP3 = table[[3]], tableP4 = table[[4]], tableP5 = table[[5]],
                tableP6 = table[[6]], tableP7 = table[[7]], tableP8 = table[[8]], tableP9 = table[[9]], tableP10 = table[[10]],
                tableP11 = table[[11]], tableP12 = table[[12]], tableP13 = table[[13]], tableP14 = table[[14]], tableP15 = table[[15]],
                tableNP1 = tableN[[1]], tableNP2 = tableN[[2]], tableNP3 = tableN[[3]], tableNP4 = tableN[[4]], tableNP5 = tableN[[5]],
                tableNP6 = tableN[[6]], tableNP7 = tableN[[7]], tableNP8 = tableN[[8]], tableNP9 = tableN[[9]], tableNP10 = tableN[[10]],
                tableNP11 = tableN[[11]], tableNP12 = tableN[[12]], tableNP13 = tableN[[13]], tableNP14 = tableN[[14]], tableNP15 = tableN[[15]],
                n = object$n, N = object$N,
                # sigma2 = object$sigma2, sigma = object$sigma2, Model = object$Model,
                # nu = object$nu, sigma2.a = object$sigma2.a, sigma.a = object$sigma2.a,
                # nu.a = object$nu.a,
                idx.lab = idx.lab,
                short.formula = object$short.formula, suStf = suStf,
                # formula1 = formula1, formula2 = formula2, formula3 = formula3, formula4 = formula4, formula5 = formula5,
                # formula6 = formula6, formula7 = formula7, formula8 = formula8, formula9 = formula9, formula10 = formula10,
                # formula11 = formula11, formula12 = formula12, formula13 = formula13, formula14 = formula14, formula15 = formula15,
                t.edf = object$msm.post.object$t.edf,
                postpost = object$msm.post.object$msm.post.object,
                state.pairs.CT.out = state.pairs.CT.out
                # CImu = CImu, mu = mu, CIsig = CIsig2, CInu = CInu, margins = object$margins,
                # l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3,
                # l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6,
                # l.sp7 = object$l.sp7, l.sp8 = object$l.sp8, X2.null = is.null(object$X2),
                # univar.gamlss = TRUE, surv.flex = object$surv.flex,
                # K1 = NULL, robust = object$robust, indx = object$fit$indx
                )
    class(res) <- "summary.fmsm"
    res

}
