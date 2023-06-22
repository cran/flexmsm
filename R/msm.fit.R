




msm.fit = function(params.0, sp, pen.matr.S.lambda, pmethod, suStf, death,
                   Q.diagnostics = FALSE, iterlimsp = iterlimsp, iterlim = iterlim,
                   sp.method = 'perf',
                   verbose, tolsp, tolsp.EFS, parallel, no_cores){

  stoprule.SP = stoprule.SP.efs = NULL
  conv.sp = Qmatr.diagnostics.list = NULL
  rp = D = L = Sl.sfTemp = St = NULL # will be not null only if efs is used
  bs.mgfit = magpp = wor.c = NULL

  # Unpack necessary variables -------------------------------------------------------------------
  data = suStf$data # HERE WE USED TO HAVE $data BUT I THINK THAT IN ALL OTHER TESTS I NORMALLY DO data AND data.long COINCIDE (not in no intercept test though)
  nstates = suStf$nstates
  start.pos.par = suStf$start.pos.par
  do.gradient = do.hessian = TRUE # not sure if to keep this or not

  pos.optparams = suStf$pos.optparams
  pos.optparams2 = suStf$pos.optparams2

  l.short.formula = suStf$l.short.formula
  whereQ = suStf$whereQ
  full.X = suStf$full.X

  Sl.sf = suStf$Sl.sf # what is the difference between Sl.sf and S.list ? (check this - it seems former is only for EFS)

  silent.trust = FALSE


  # Fixing S.list
  S.list = suStf$S.list
  # Keep an eye on whether this fix is OK - this is just to remove NULL, may not be necessary any more
  S.list.temp = list()
  ii.fix = 1
  for(i.fix in 1:length(S.list)){
    if(!is.null(S.list[[i.fix]])){
      S.list.temp[[ii.fix]] = S.list[[i.fix]]
      ii.fix = ii.fix + 1
    }
  }
  S.list = S.list.temp
  rm(S.list.temp)

  start.pos.par.only.smooth = suStf$start.pos.par.only.smooth
  start.pos.par.only.smooth.FPC = suStf$start.pos.par.only.smooth.FPC
  start.pos.par.detailed = suStf$start.pos.par.detailed

  MM = list(start.pos.par = start.pos.par,
            pos.optparams = pos.optparams,
            pos.optparams2 = pos.optparams2,
            l.short.formula = l.short.formula,
            whereQ = whereQ,
            nstates = nstates,
            l.params = length(params.0),
            cens.state = suStf$cens.state)


  # ******************** #
  # Penalty matrix setup #
  # ******************** #
  # Setup full penalty matrix to be used for penalized likelihood estimation
  pen.matr.S.lambda = penalty.setup(sp = sp, suStf = suStf)
  pen.matr.S.lambda = pen.matr.S.lambda[1:max(pos.optparams2) == pos.optparams2, 1:max(pos.optparams2) == pos.optparams2]
  # ----------------------------------------------------------------------------------------------


  # ***** OPTIMIZATION STARTS HERE *****

  if(Q.diagnostics) Qmatr.diagnostics.list = list()

  ii = 1
  fit.all = list()

  fit <- try(trust::trust(LikGradHess.general, params.0, rinit = 1, # the following are all params needed by the function ********
             data = data, full.X = full.X, MM = MM, pen.matr.S.lambda = pen.matr.S.lambda, # *****
              aggregated.provided = FALSE, do.gradient = TRUE, do.hessian = TRUE,
             pmethod = pmethod, death = death,
             Qmatr.diagnostics.list = Qmatr.diagnostics.list,
             parallel = parallel, no_cores = no_cores,
             rmax = 100, parscale = rep(1, length(params.0)), blather = TRUE, iterlim = iterlim), silent = silent.trust)

  if(!is.null(attr(fit$error, 'class'))) warning('Unable to carry out call n.', ii, 'to trust region algorithm. Look into partial results for diagnostics.')



  fit.all[[ii]] = fit
  ii = ii + 1

  sp.all = c(sp)


  if(length(sp) > 0){ # i.e. only do smoothing parameters estimation if there are smooths! (change/check condition)

    # ********************************************* #
    # MAGIC BASED SMOOTHING PARAMETER ESTIMATION ####
    # ********************************************* #
    if(sp.method == 'perf'){

      stoprule.SP = 1
      tolsp = tolsp
      iter.sp = 0
      iter.inner = 0

      while (stoprule.SP > tolsp) {

        # Saving starting parameters (lambda and beta) and initial -log-likelihood
        fito <- fit$l
        o.ests <- c(fit$argument)

        spo <- sp

        # Obtaining necessay quantities for magic function (i.e. lambda optimization)
        wor.c <- GJRM::working.comp(fit)

        bs.mgfit <- mgcv::magic(y = wor.c$Z, X = wor.c$X,
                          sp = sp, S = S.list, off = suStf$start.pos.par.only.smooth.FPC.constr,
                          # rank = qu.mag$rank,
                          gcv = FALSE)
        sp <- bs.mgfit$sp

        iter.sp <- iter.sp + 1
        names(sp) <- names(spo)

        # Setup full penalty matrix to be used for penalized likelihood estimation
        pen.matr.S.lambda = penalty.setup(sp = sp, suStf = suStf)
        pen.matr.S.lambda = pen.matr.S.lambda[1:max(pos.optparams2) == pos.optparams2, 1:max(pos.optparams2) == pos.optparams2]


        # ***** With Q computation put internally *****
        if(Q.diagnostics) Qmatr.diagnostics.list = fit$Qmatr.diagnostics.list

        fit <- try(trust::trust(LikGradHess.general, o.ests, rinit = 1, # the following are all params needed by the function ********
                     data = data, full.X = full.X, MM = MM, pen.matr.S.lambda = pen.matr.S.lambda, # *****
                      aggregated.provided = FALSE, do.gradient = TRUE, do.hessian = TRUE,
                     pmethod = pmethod, death = death,
                     Qmatr.diagnostics.list = Qmatr.diagnostics.list,
                     parallel = parallel, no_cores = no_cores, # **************************************
                     rmax = 100, parscale = rep(1, length(params.0)), blather = TRUE, iterlim = iterlim), silent = silent.trust)

        if(!is.null(attr(fit$error, 'class'))) warning('Unable to carry out call n.', ii, 'to trust region algorithm. Look into partial results for diagnostics.')

        fit.all[[ii]] = fit
        ii = ii + 1

        iter.inner <- iter.inner + fit$iterations
        sp.all = c(sp.all, sp)

        if (iter.sp > iterlimsp) {
          conv.sp <- FALSE
          warning('Maximum number of iterations for smoothing parameter estimation reached.')
          break
        }

        stoprule.SP <- abs(fit$l - fito)/(0.1 + abs(fit$l))
        # if(verbose) print(stoprule.SP)

      }


      magpp <- mgcv::magic.post.proc(wor.c$X, bs.mgfit) # SHOULD wor.c and bs.mgfit ALSO BE UPDATED WITH THE NEW FIT ??? GIAMPIERO CODE DOES NOT

    }


    if(sp.method == 'efs'){ # EFS SMOOTHING ####

      # ******************************************* #
      # EFS BASED SMOOTHING PARAMETER ESTIMATION ####
      # ******************************************* #


      qu.mag = list(Ss = S.list,
                    off =  suStf$start.pos.par.only.smooth.FPC.constr)


      LDfun <- function(Hp, eigen.fix) {
        rank <- dim(Hp)[1]
        D <- diag(Hp)
        if (sum(!is.finite(D)) > 0)
          stop("non finite values in Hessian")
        if (min(D) < 0) {
          Dthresh <- max(D) * sqrt(.Machine$double.eps)
          if (-min(D) < Dthresh) {
            indefinite <- FALSE
            D[D < Dthresh] <- Dthresh
          }
          else indefinite <- TRUE
        }
        else indefinite <- FALSE
        if (indefinite) {
          if (eigen.fix) {
            eh <- eigen(Hp, symmetric = TRUE)
            ev <- abs(eh$values)
            Hp <- eh$vectors %*% (ev * t(eh$vectors))
          }
          else {
            Ib <- diag(rank) * abs(min(D))
            Ip <- diag(rank) * abs(max(D) * .Machine$double.eps^0.5)
            Hp <- Hp + Ip + Ib
          }
          D <- rep(1, ncol(Hp))
          indefinite <- TRUE
        }
        else {
          D <- D^-0.5
          Hp <- D * t(D * Hp)
          Ip <- diag(rank) * .Machine$double.eps^0.5
        }
        L <- suppressWarnings(chol(Hp, pivot = TRUE))
        mult <- 1
        while (attr(L, "rank") < rank) {
          if (eigen.fix) {
            eh <- eigen(Hp, symmetric = TRUE)
            ev <- eh$values
            thresh <- max(min(ev[ev > 0]), max(ev) *
                            1e-06) * mult
            mult <- mult * 10
            ev[ev < thresh] <- thresh
            Hp <- eh$vectors %*% (ev * t(eh$vectors))
            L <- suppressWarnings(chol(Hp, pivot = TRUE))
          }
          else {
            L <- suppressWarnings(chol(Hp + Ip, pivot = TRUE))
            Ip <- Ip * 100
          }
          indefinite <- TRUE
        }
        list(L = L, D = D)
      }



      stoprule.SP <- 1
      conv.sp <- TRUE
      iter.inner <- iter.sp <- 0
      controlEFS <- list(efs.lspmax = 15, eps = 1e-7, # this is the one used when iter == 1, can leave like this
                         tol = 1e-06, tiny = .Machine$double.eps^0.5,
                         efs.tol = tolsp.EFS)
      score.hist <- rep(0, 200)
      mult <- 1
      lsp <- log(sp)
      gamma <- 1
      Mp <- -1
      eigen.fix <- FALSE

      Sl.termMult <- getFromNamespace("Sl.termMult", "mgcv")
      ldetS <- getFromNamespace("ldetS", "mgcv")

      Mp <- ncol(totalPenaltySpace(qu.mag$Ss, NULL, qu.mag$off,
                                   length(fit$argument))$Z)


      # Fixing Sl.sf
      Sl.sf.temp = list()
      jj = 1
      smooth.mapping = match(start.pos.par.detailed[-length(start.pos.par.detailed)], start.pos.par.only.smooth[-length(start.pos.par.only.smooth)])
      smooth.length = diff(start.pos.par.detailed)[!is.na(smooth.mapping)]
      for(i in 1:length(Sl.sf)){

        for(i.intern in 1:length(Sl.sf[[i]])) {
          Sl.sf.temp[[jj]] = Sl.sf[[i]][[i.intern]]

          Sl.sf.temp[[jj]]$start = suStf$start.pos.par.only.smooth.constr[[jj]]
          Sl.sf.temp[[jj]]$stop = suStf$start.pos.par.only.smooth.constr[[jj]] + smooth.length[jj] - 1 # IS THIS OK ALSO FOR MULTIPLE SMOOTHS per transition ???

          jj = jj + 1
        }
      }

      Sl.sf = Sl.sf.temp
      rm(Sl.sf.temp)

      # attr(Sl.sf, "lambda") <- attr(Sl.sf, "lambda") # not sure if these two are needed too - leave for now
      # attr(Sl.sf, "E") <- attr(Sl.sf, "E")
      attr(Sl.sf, 'cholesky') = FALSE


      Sl.sfTemp <- Sl.sf
      for (i in 1:length(Sl.sfTemp)){
        Sl.sfTemp[[i]]$D <- solve(Sl.sfTemp[[i]]$D)
      }
      # ***********************************************


      for (iter in 1:200) {
        o.ests <- c(fit$argument)
        rp <- ldetS(Sl.sf, rho = lsp, fixed = rep(FALSE, length(lsp)), np = length(fit$argument), root = TRUE)
        o.estsStar <- Sl.initial.repara(Sl.sfTemp, o.ests, inverse = TRUE)

        Vb <- Sl.initial.repara(Sl.sfTemp, GJRM::PDef(fit$hessian)$res.inv, inverse = TRUE)

        LD <- LDfun(GJRM::PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D
        ipiv <- piv <- attr(L, "pivot")
        p <- length(piv)
        ipiv[piv] <- 1:p
        Vb <- crossprod(forwardsolve(t(L), diag(D, nrow = p)[piv, , drop = FALSE])[ipiv, , drop = FALSE])
        Vb <- Sl.repara(rp$rp, Vb, inverse = TRUE)
        SVb <- Sl.termMult(Sl.sf, Vb)
        trVS <- rep(0, length(SVb))

        for (i in 1:length(SVb)) {
          ind <- attr(SVb[[i]], "ind")
          trVS[i] <- sum(diag(SVb[[i]][, ind]))
        }

        start <- Sl.repara(rp$rp, o.estsStar)
        Sb <- Sl.termMult(Sl.sf, start, full = TRUE)
        bSb <- rep(0, length(Sb))
        for (i in 1:length(Sb)) bSb[i] <- sum(start * Sb[[i]])

        S1 <- rp$ldet1
        a <- pmax(controlEFS$tiny, S1 * exp(-lsp) - trVS)
        r <- a/pmax(controlEFS$tiny, bSb)
        r[a == 0 & bSb == 0] <- 1
        r[!is.finite(r)] <- 1e+06
        lsp1 <- pmin(lsp + log(r) * mult, controlEFS$efs.lspmax)
        max.step <- max(abs(lsp1 - lsp))
        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        old.reml <- -as.numeric((-fit$l - drop(t(fit$argument) %*% fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)

        sp1 <- exp(lsp1)

        names(sp1) <- names(lsp1)


        # Setup full penalty matrix to be used for penalized likelihood estimation
        pen.matr.S.lambda = penalty.setup(sp = sp1, suStf = suStf)
        pen.matr.S.lambda = pen.matr.S.lambda[1:max(pos.optparams2) == pos.optparams2, 1:max(pos.optparams2) == pos.optparams2]


        if(Q.diagnostics) Qmatr.diagnostics.list = fit$Qmatr.diagnostics.list
        fit <- try(trust::trust(LikGradHess.general, o.ests, rinit = 1, # the following are all params needed by the function ********
                     data = data, full.X = full.X, MM = MM, pen.matr.S.lambda = pen.matr.S.lambda, # *****
                      aggregated.provided = FALSE, do.gradient = TRUE, do.hessian = TRUE,
                     pmethod = pmethod, death = death,
                     Qmatr.diagnostics.list = Qmatr.diagnostics.list,
                     parallel = parallel, no_cores = no_cores, # *****
                     rmax = 100, parscale = rep(1, length(params.0)), blather = TRUE, iterlim = iterlim), silent = silent.trust)


        if(!is.null(attr(fit$error, 'class'))) warning('Unable to carry out call n.', ii, 'to trust region algorithm. Look into partial results for diagnostics.')

        iter.inner <- iter.inner + fit$iterations
        rp <- ldetS(Sl.sf, rho = lsp1, fixed = rep(FALSE, length(lsp1)), np = length(fit$argument), root = TRUE)
        Vb <- Sl.initial.repara(Sl.sfTemp, GJRM::PDef(fit$hessian)$res.inv, inverse = TRUE)
        LD <- LDfun(GJRM::PDef(Vb)$res.inv, eigen.fix)
        L <- LD$L
        D <- LD$D
        ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
        fit$REML <- -as.numeric((-fit$l - drop(t(fit$argument) %*% fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)

        if (fit$REML <= old.reml) { # ***** NEED TO RUN FROM THIS LINE AND ON *****
          # print('I made the REML decrease')
          if (max.step < 0.05) {
            lsp2 <- pmin(lsp + log(r) * mult * 2, 12)
            sp2 <- exp(lsp2)
            names(sp2) <- names(lsp2)

            # Setup full penalty matrix to be used for penalized likelihood estimation
            pen.matr.S.lambda = penalty.setup(sp = sp2, suStf = suStf)
            pen.matr.S.lambda = pen.matr.S.lambda[1:max(pos.optparams2) == pos.optparams2, 1:max(pos.optparams2) == pos.optparams2]

            if(Q.diagnostics) Qmatr.diagnostics.list = fit$Qmatr.diagnostics.list
            fit2 <- try(trust::trust(LikGradHess.general, o.ests, rinit = 1, # the following are all params needed by the function ********
                          data = data, full.X = full.X, MM = MM, pen.matr.S.lambda = pen.matr.S.lambda, # *****
                           aggregated.provided = FALSE, do.gradient = TRUE, do.hessian = TRUE,
                          pmethod = pmethod, death = death,
                          Qmatr.diagnostics.list = Qmatr.diagnostics.list,
                          parallel = parallel, no_cores = no_cores, # *****
                          rmax = 100, parscale = rep(1, length(params.0)), blather = TRUE, iterlim = iterlim), silent = silent.trust)


            if(!is.null(attr(fit2$error, 'class'))) warning('Unable to carry out call n.', ii, 'to trust region algorithm. Look into partial results for diagnostics.')


            iter.inner <- iter.inner + fit2$iterations
            rp <- ldetS(Sl.sf, rho = lsp2, fixed = rep(FALSE, length(lsp2)), np = length(fit$argument), root = TRUE)
            Vb <- Sl.initial.repara(Sl.sfTemp, GJRM::PDef(fit2$hessian)$res.inv,
                                    inverse = TRUE)
            LD <- LDfun(GJRM::PDef(Vb)$res.inv, eigen.fix)
            L <- LD$L
            D <- LD$D
            ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
            fit2$REML <- -as.numeric((-fit2$l - drop(t(fit2$argument) %*% fit2$S.h %*% fit2$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)

            if (fit2$REML < fit$REML) {
              fit <- fit2
              lsp <- lsp2
              # print('I am using lsp2 because REML decreased, max step was < 0.05 and sp2 is better so increased the mult')
              mult <- mult * 2
            } else {
              lsp <- lsp1
              # print('I am using lsp1 because REML decreased, max step was < 0.05 but sp2 is not better')
            }
          } else {
            lsp <- lsp1
            # print('I am using lsp1 because REML decreased and the max step was >= 0.05')
          }
        } else {
          jj = 0
          while (fit$REML > old.reml && mult > 1) {
            mult <- mult/2

            jj = jj + 1 # just for debug

            lsp1 <- pmin(lsp + log(r) * mult, controlEFS$efs.lspmax)
            sp1 <- exp(lsp1)

            names(sp1) <- names(lsp1)

            # Setup full penalty matrix to be used for penalized likelihood estimation
            pen.matr.S.lambda = penalty.setup(sp = sp1, suStf = suStf)
            pen.matr.S.lambda = pen.matr.S.lambda[1:max(pos.optparams2) == pos.optparams2, 1:max(pos.optparams2) == pos.optparams2]


            if(Q.diagnostics) Qmatr.diagnostics.list = fit$Qmatr.diagnostics.list

            fit <- try(trust::trust(LikGradHess.general, o.ests, rinit = 1, # the following are all params needed by the function ********
                         data = data, full.X = full.X, MM = MM, pen.matr.S.lambda = pen.matr.S.lambda, # *****
                          aggregated.provided = FALSE, do.gradient = TRUE, do.hessian = TRUE,
                         pmethod = pmethod, death = death,
                         Qmatr.diagnostics.list = Qmatr.diagnostics.list,
                         parallel = parallel, no_cores = no_cores, # *****
                         rmax = 100, parscale = rep(1, length(params.0)), blather = TRUE, iterlim = iterlim), silent = silent.trust)


            if(!is.null(attr(fit$error, 'class'))) warning('Unable to carry out call n.', ii, 'to trust region algorithm. Look into partial results for diagnostics.')

            iter.inner <- iter.inner + fit$iterations
            rp <- ldetS(Sl.sf, rho = lsp1, fixed = rep(FALSE, length(lsp1)), np = length(fit$argument), root = TRUE)
            Vb <- Sl.initial.repara(Sl.sfTemp, GJRM::PDef(fit$hessian)$res.inv, inverse = TRUE)
            LD <- LDfun(GJRM::PDef(Vb)$res.inv, eigen.fix)
            L <- LD$L
            D <- LD$D
            ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
            fit$REML <- -as.numeric((-fit$l - drop(t(fit$argument) %*% fit$S.h %*% fit$argument)/2)/gamma + rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)
          }
          lsp <- lsp1
          # if(verbose) print(paste('I am using alternative lsp1 because REML did not decrease - I shrunk the mult', jj, 'times.'))

          if (mult < 1) mult <- 1
        }
        score.hist[iter] <- fit$REML
        if (iter > 3 && max.step < 0.05 && max(abs(diff(score.hist[(iter - 3):iter]))) < controlEFS$efs.tol){

          break
        }
        if (iter == 1){
          old.ll <- fit$l
        } else {
          stoprule.SP.efs = abs(old.ll - fit$l) / (100* controlEFS$eps * abs(fit$l))
          if (abs(old.ll - fit$l) < 100 * controlEFS$eps * abs(fit$l)) break
          old.ll <- fit$l
        }

        # if(verbose) print(exp(lsp))

        sp.all = c(sp.all, exp(lsp))
        fit.all[[ii]] = fit
        ii = ii + 1

      }

      sp <- exp(lsp)
      iter.sp <- iter

      if (iter > 200){
        conv.sp <- FALSE
      } else {
        conv.sp <- TRUE
      }

      St <- crossprod(rp$E)
    }


  } else {
  wor.c <- GJRM::working.comp(fit)
  bs.mgfit <- mgcv::magic(y = wor.c$Z, X = wor.c$X,
                    sp = numeric(0), S = list(),
                    off = numeric(0))
  magpp <- mgcv::magic.post.proc(wor.c$X, bs.mgfit)
}



  sp.all.matr = matrix(sp.all, nrow = length(sp))


  list(fit = fit,
       fit.all = fit.all, sp.all = sp.all, sp.all.matr = sp.all.matr, # later on will only keep the matrix version, but need to make sure the definition is OK
       iter.if = fit.all[[1]]$iterations,
       sp.method = sp.method, conv.sp = conv.sp, iter.sp = iter.sp,
       iter.inner = iter.inner, bs.mgfit = bs.mgfit, wor.c = wor.c,
       sp = sp, magpp = magpp, Sl.sf = Sl.sf,
       Sl.sfTemp = Sl.sfTemp, # from efs, used later for post fit quantities computation
       rp = rp$rp, D = D, L = L, St = St,
       Qmatr.diagnostics.list = Qmatr.diagnostics.list,
       stoprule.SP = stoprule.SP, stoprule.SP.efs = stoprule.SP.efs)




}
