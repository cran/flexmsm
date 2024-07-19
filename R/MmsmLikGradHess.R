


MmsmLikGradHess = function(params,
                           formula1 = formula1, data1 = data1, l.params1 = length(params1), # --- arguments needed for the objfun from here ---
                           id1 = id1, state1 = state1, spP1 = spP1,
                           formula2 = formula2, data2 = data2, l.params2 = length(params2),
                           id2 = id2, state2 = state2,  spP2 = spP2,
                           pmethod = pmethod,
                           Q.diagnostics = Q.diagnostics, iterlimsp = iterlimsp,
                           iterlim = iterlim,
                           sp.method = sp.method,
                           verbose = verbose, tolsp = tolsp, tolsp.EFS = tolsp.EFS,
                           parallel = parallel, no_cores = no_cores){


  if(length(params) != (l.params1 + l.params2 + 1)) stop(paste('l.params1 and/or l.params2 are incorrect:',
                                                               length(params), 'parameters',
                                                               'were provided but (l.params1 + l.params2 + 1) =',
                                                               (l.params1 + l.params2 + 1) ))

  params1 = params[1:l.params1]
  params2 = params[(l.params1+1):(l.params1+l.params2)]
  phi = exp(params[l.params1+l.params2+1])


  # ********************** #
  # SETUP for PROCESS 1 ####
  # ********************** #

  proc1 <- try(fmsm(formula = formula1, data = data1,
                    params.0 = params1, sp.0 = spP1,
                    id = id1, state = state1, death = FALSE,
                    fit = FALSE, justComp = c('lik', 'grad', 'hess'),
                    pmethod = pmethod,
                    Q.diagnostics = Q.diagnostics, iterlimsp = iterlimsp,
                    iterlim = iterlim,
                    sp.method = sp.method,
                    verbose = verbose, tolsp = tolsp, tolsp.EFS = tolsp.EFS,
                    parallel = parallel, no_cores = no_cores),
               silent = TRUE)


  # *********************** #
  # SETUP for PROCESS 2  ####
  # *********************** #

  proc2 <- try(fmsm(formula = formula2, data = data2,
                    params.0 = params2, sp.0 = spP2,
                    id = id2, state = state2, death = FALSE,
                    fit = FALSE, justComp = c('lik', 'grad', 'hess'),
                    pmethod = pmethod,
                    Q.diagnostics = Q.diagnostics, iterlimsp = iterlimsp,
                    iterlim = iterlim,
                    sp.method = sp.method,
                    verbose = verbose, tolsp = tolsp, tolsp.EFS = tolsp.EFS,
                    parallel = parallel, no_cores = no_cores),
               silent = TRUE)


  # ******************************************** #
  # FULL LOG-LIKELIHOOD, GRADIENT AND HESSIAN ####
  # ******************************************** #

  nstates = proc1$suStf$nstates
  copContr = 0
  copContrvec = c()

  # initialise gradient quantities
  dcopContr.d1 = dcopContr.d2 = dcopContr.dphi = 0

  # initialise hessian quantities
  d2copContr.d1d1 = d2copContr.d1d2 = d2copContr.d1dphi = 0
  d2copContr.d2d2 = d2copContr.d2dphi = d2copContr.dphidphi = 0

  # for debugging purpose
  Vecd2copContr.d1d1 = Vecd2copContr.d1d2 = Vecd2copContr.d1dphi = c()
  Vecd2copContr.d2d2 = Vecd2copContr.d2dphi = Vecd2copContr.dphidphi = c()

  l.params1 = length(proc1$suStf$params)
  l.params2 = length(proc2$suStf$params)

  jj.marg1 = jj.marg2 = 0

  for(inds in 1:length(proc1$singleComp$P.CM.contr)
  ){ # anyway there will also be the same number of *individuals* in both processes (by structure)

    # print(inds) # needed for debugging only

    # EXTRACTION OF THE MARGINAL TERMS - PROCESS 1 ####
    Ps1 = diag(nrow = proc1$suStf$nstates)
    dPs1 = array(0, dim = dim(proc1$singleComp$dP.hist)[c(1,2,3)])
    d2Ps1 = array(0, dim = dim(proc1$singleComp$d2P.hist)[c(1,2,3,4)])

    if(proc1$singleComp$ind.CM$ni[inds] > 1){

      for(i.m1 in 1:(proc1$singleComp$ind.CM$ni[inds]-1)){ # this is ni-1 in case we have the R time (added later, separately)

        # P based term
        Ps1 = Ps1 %*% proc1$singleComp$P.hist[,,jj.marg1 + i.m1]

        # dP based term
        PprodsBefore = PprodsAfter = diag(nstates)

        if(i.m1 > 1) for(j in 1:(i.m1-1)) PprodsBefore = PprodsBefore %*% proc1$singleComp$P.hist[,,jj.marg1 + j]
        if(i.m1 < (proc1$singleComp$ind.CM$ni[inds]-1)) for(j in (i.m1+1):(proc1$singleComp$ind.CM$ni[inds]-1)) PprodsAfter = PprodsAfter %*% proc1$singleComp$P.hist[,,jj.marg1 + j]

        for(lpar in 1:l.params1){
          dPs1[,,lpar] = dPs1[,,lpar] + PprodsBefore %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1 + i.m1] %*% PprodsAfter
          for(l2par in 1:l.params1){
            if(l2par >= lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,lpar,l2par,jj.marg1 + i.m1] # upper triangle: OK, we have these in d2P.hist
            if(l2par < lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,l2par,lpar,jj.marg1 + i.m1]  # lower triangle: take elements from upper triangle and flip them
            d2Ps1[,,lpar,l2par] = d2Ps1[,,lpar,l2par] + PprodsBefore %*% d2Pcmpnt %*% PprodsAfter
          }
        }

        # d2P based term
        if(i.m1 < (proc1$singleComp$ind.CM$ni[inds]-1)){
          for(i2.m1 in (i.m1+1):(proc1$singleComp$ind.CM$ni[inds]-1)){

            PprodsBefore2 = PprodsAfter2 = diag(nstates)

            if(i2.m1 > (i.m1+1)) for(j in (i.m1+1):(i2.m1-1)) PprodsBefore2 = PprodsBefore2 %*% proc1$singleComp$P.hist[,,jj.marg1 + j]
            if(i2.m1 < (proc1$singleComp$ind.CM$ni[inds]-1)) for(j in (i2.m1+1):(proc1$singleComp$ind.CM$ni[inds]-1)) PprodsAfter2 = PprodsAfter2 %*% proc1$singleComp$P.hist[,,jj.marg1 + j]

            for(lpar in 1:l.params1){
              beforeX = PprodsBefore %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1 + i.m1]
              afterXsw = PprodsBefore2 %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1 + i2.m1] %*% PprodsAfter2
              for(l2par in 1:l.params1){
                afterX = PprodsBefore2 %*% proc1$singleComp$dP.hist[,,l2par,jj.marg1 + i2.m1] %*% PprodsAfter2
                beforeXsw = PprodsBefore %*% proc1$singleComp$dP.hist[,,l2par,jj.marg1 + i.m1]

                d2Ps1[,,lpar,l2par] = d2Ps1[,,lpar,l2par] + beforeX %*% afterX + beforeXsw %*% afterXsw
              }
            }
          }
        }

        # print(d2Ps1[1,3,,])  # for debugging purposes
      }
    }
    jj.marg1 = jj.marg1 + proc1$singleComp$ind.CM$ni[inds]


    # EXTRACTION OF THE MARGINAL TERMS - PROCESS 2 ####
    Ps2 = diag(nrow = proc2$suStf$nstates)
    dPs2 = array(0, dim = dim(proc2$singleComp$dP.hist)[c(1,2,3)])
    d2Ps2 = array(0, dim = dim(proc2$singleComp$d2P.hist)[c(1,2,3,4)])

    if(proc2$singleComp$ind.CM$ni[inds] > 1){

      for(i.m2 in 1:(proc2$singleComp$ind.CM$ni[inds]-1)){ # this is ni-1 in case we have the R time (added later, separately)

        # P based term
        Ps2 = Ps2 %*% proc2$singleComp$P.hist[,,jj.marg2 + i.m2]

        # dP based term
        PprodsBefore = PprodsAfter = diag(nstates)

        if(i.m2 > 1) for(j in 1:(i.m2-1)) PprodsBefore = PprodsBefore %*% proc2$singleComp$P.hist[,,jj.marg2 + j]
        if(i.m2 < (proc2$singleComp$ind.CM$ni[inds]-1)) for(j in (i.m2+1):(proc2$singleComp$ind.CM$ni[inds]-1)) PprodsAfter = PprodsAfter %*% proc2$singleComp$P.hist[,,jj.marg2 + j]


        for(lpar in 1:l.params2){
          dPs2[,,lpar] = dPs2[,,lpar] + PprodsBefore %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2 + i.m2] %*% PprodsAfter
          for(l2par in 1:l.params2){
            if(l2par >= lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,lpar,l2par,jj.marg2 + i.m2] # upper triangle: OK, we have these in d2P.hist
            if(l2par < lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,l2par,lpar,jj.marg2 + i.m2]  # lower triangle: take elements from upper triangle and flip them
            d2Ps2[,,lpar,l2par] = d2Ps2[,,lpar,l2par] + PprodsBefore %*% d2Pcmpnt %*% PprodsAfter
          }
        }

        # d2P based term
        if(i.m2 < (proc2$singleComp$ind.CM$ni[inds]-1)){
          for(i2.m2 in (i.m2+1):(proc2$singleComp$ind.CM$ni[inds]-1)){

            PprodsBefore2 = PprodsAfter2 = diag(nstates)

            if(i2.m2 > (i.m2+1)) for(j in (i.m2+1):(i2.m2-1)) PprodsBefore2 = PprodsBefore2 %*% proc2$singleComp$P.hist[,,jj.marg2 + j]
            if(i2.m2 < (proc2$singleComp$ind.CM$ni[inds]-1)) for(j in (i2.m2+1):(proc2$singleComp$ind.CM$ni[inds]-1)) PprodsAfter2 = PprodsAfter2 %*% proc2$singleComp$P.hist[,,jj.marg2 + j]

            for(lpar in 1:l.params2){
              beforeX = PprodsBefore %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2 + i.m2]
              afterXsw = PprodsBefore2 %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2 + i2.m2] %*% PprodsAfter2
              for(l2par in 1:l.params2){
                afterX = PprodsBefore2 %*% proc2$singleComp$dP.hist[,,l2par,jj.marg2 + i2.m2] %*% PprodsAfter2
                beforeXsw = PprodsBefore %*% proc2$singleComp$dP.hist[,,l2par,jj.marg2 + i.m2]

                d2Ps2[,,lpar,l2par] = d2Ps2[,,lpar,l2par] + beforeX %*% afterX + beforeXsw %*% afterXsw
              }
            }
          }
        }
      }
    }
    jj.marg2 = jj.marg2 + proc2$singleComp$ind.CM$ni[inds]



    # COMPOSITE MODEL LIKELIHOOD, GRADIENT AND HESSIAN ####

    if(proc1$singleComp$ind.CM$obsT[inds] == 2 & proc2$singleComp$ind.CM$obsT[inds] == 2){
      m1L = 1 - Ps1[1,nstates]
      m2L = 1 - Ps2[1,nstates]

      m1R = 1 - (Ps1 %*% proc1$singleComp$P.hist[,,jj.marg1])[1,nstates]
      m2R = 1 - (Ps2 %*% proc2$singleComp$P.hist[,,jj.marg2])[1,nstates]

      # LIKELIHOOD
      copContri = copFun(m1R, m2R, phi, type = 'C0') + copFun(m1L, m2L, phi, type = 'C0') - copFun(m1R, m2L, phi, type = 'C0') - copFun(m1L, m2R, phi, type = 'C0') # cross contribs and two right times
      copContr = copContr + log(copContri)
      copContrvec = c(copContrvec, copContri)

      # GRADIENT (we will have three vectors here: dtheta1, dtheta2, dphi)
      # dtheta1
      dm1L = -dPs1[1,nstates,]

      dm1R = array(0, dim(proc1$singleComp$dP.hist)[3])
      for(lpar in 1:l.params1) dm1R[lpar] = -(dPs1[,,lpar] %*% proc1$singleComp$P.hist[,,jj.marg1] + Ps1 %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1])[1,nstates]

      dDM.d1 = ( dcopFun.d1(m1R, m2R, phi, type = 'C0') * dm1R + dcopFun.d1(m1L, m2L, phi, type = 'C0') * dm1L
                 - dcopFun.d1(m1R, m2L, phi, type = 'C0') * dm1R - dcopFun.d1(m1L, m2R, phi, type = 'C0') * dm1L)
      dcopContr.d1 = dcopContr.d1 + copContri^-1 * dDM.d1

      # dtheta2
      dm2L = -dPs2[1,nstates,]

      dm2R = array(0, dim(proc2$singleComp$dP.hist)[3])
      for(lpar in 1:l.params2) dm2R[lpar] = -(dPs2[,,lpar] %*% proc2$singleComp$P.hist[,,jj.marg2] + Ps2 %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2])[1,nstates]

      dDM.d2 = ( dcopFun.d2(m1R, m2R, phi, type = 'C0') * dm2R + dcopFun.d2(m1L, m2L, phi, type = 'C0') * dm2L
                 - dcopFun.d2(m1R, m2L, phi, type = 'C0') * dm2L - dcopFun.d2(m1L, m2R, phi, type = 'C0') * dm2R)
      dcopContr.d2 = dcopContr.d2 + copContri^-1 * dDM.d2

      #dphi
      dDM.dphi = (dcopFun.dphi(m1R, m2R, phi, type = 'C0') + dcopFun.dphi(m1L, m2L, phi, type = 'C0')
                  - dcopFun.dphi(m1R, m2L, phi, type = 'C0') - dcopFun.dphi(m1L, m2R, phi, type = 'C0')) * phi
      dcopContr.dphi = dcopContr.dphi + copContri^-1 * dDM.dphi

      # HESSIAN

      # computation of d2P terms
      # (process 1)
      d2m1L = -d2Ps1[1,nstates,,]

      d2m1R = array(0, dim(proc1$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params1){
        for(l2par in 1:l.params1){
          if(l2par >= lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,lpar,l2par,jj.marg1] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,l2par,lpar,jj.marg1]  # lower triangle: take elements from upper triangle and flip them
          d2m1R[lpar, l2par] = -(d2Ps1[,,lpar,l2par] %*% proc1$singleComp$P.hist[,,jj.marg1] + dPs1[,,l2par] %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1]
                                 + dPs1[,,lpar] %*% proc1$singleComp$dP.hist[,,l2par,jj.marg1] + Ps1 %*% d2Pcmpnt)[1,nstates]
        }
      }


      # (process 2)
      d2m2L = -d2Ps2[1,nstates,,]

      d2m2R = array(0, dim(proc2$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params2){
        for(l2par in 1:l.params2){
          if(l2par >= lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,lpar,l2par,jj.marg2] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,l2par,lpar,jj.marg2]  # lower triangle: take elements from upper triangle and flip them
          d2m2R[lpar, l2par] = -(d2Ps2[,,lpar,l2par] %*% proc2$singleComp$P.hist[,,jj.marg2] + dPs2[,,l2par] %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2]
                                 + dPs2[,,lpar] %*% proc2$singleComp$dP.hist[,,l2par,jj.marg2] + Ps2 %*% d2Pcmpnt)[1,nstates]
        }
      }

      # computation of DM derivatives and corresponding copContr
      # dtheta1.dtheta1
      d2DM.d1d1 = ( d2copFun.d1d1(m1R, m2R, phi, type = 'C0') * (dm1R %*% t(dm1R)) + dcopFun.d1(m1R, m2R, phi, type = 'C0') * d2m1R
                    + d2copFun.d1d1(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm1L)) + dcopFun.d1(m1L, m2L, phi, type = 'C0') * d2m1L
                    - d2copFun.d1d1(m1R, m2L, phi, type = 'C0') * (dm1R %*% t(dm1R)) - dcopFun.d1(m1R, m2L, phi, type = 'C0') * d2m1R
                    - d2copFun.d1d1(m1L, m2R, phi, type = 'C0') * (dm1L %*% t(dm1L)) - dcopFun.d1(m1L, m2R, phi, type = 'C0') * d2m1L
      )
      d2copContr.d1d1 = d2copContr.d1d1 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d1))  + copContri^-1 * d2DM.d1d1


      # dtheta1.dtheta2
      d2DM.d1d2 = ( d2copFun.d1d2(m1R, m2R, phi, type = 'C0') * (dm1R %*% t(dm2R))
                    + d2copFun.d1d2(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm2L))
                    - d2copFun.d1d2(m1R, m2L, phi, type = 'C0') * (dm1R %*% t(dm2L))
                    - d2copFun.d1d2(m1L, m2R, phi, type = 'C0') * (dm1L %*% t(dm2R))
      )
      d2copContr.d1d2 = d2copContr.d1d2 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d2)) + copContri^-1 * d2DM.d1d2


      # dtheta1.dphi
      d2DM.d1dphi = ( d2copFun.d1dphi(m1R, m2R, phi, type = 'C0') * dm1R
                      + d2copFun.d1dphi(m1L, m2L, phi, type = 'C0') * dm1L
                      - d2copFun.d1dphi(m1R, m2L, phi, type = 'C0') * dm1R
                      - d2copFun.d1dphi(m1L, m2R, phi, type = 'C0') * dm1L
      ) * phi
      d2DM.d1dphi = matrix(d2DM.d1dphi, nrow = length(d2DM.d1dphi))
      d2copContr.d1dphi = d2copContr.d1dphi + (-copContri^-2) * (dDM.d1 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d1dphi


      # dtheta2.dtheta2
      d2DM.d2d2 = ( d2copFun.d2d2(m1R, m2R, phi, type = 'C0') * (dm2R %*% t(dm2R)) + dcopFun.d2(m1R, m2R, phi, type = 'C0') * d2m2R
                    + d2copFun.d2d2(m1L, m2L, phi, type = 'C0') * (dm2L %*% t(dm2L)) + dcopFun.d2(m1L, m2L, phi, type = 'C0') * d2m2L
                    - d2copFun.d2d2(m1R, m2L, phi, type = 'C0') * (dm2L %*% t(dm2L)) - dcopFun.d2(m1R, m2L, phi, type = 'C0') * d2m2L
                    - d2copFun.d2d2(m1L, m2R, phi, type = 'C0') * (dm2R %*% t(dm2R)) - dcopFun.d2(m1L, m2R, phi, type = 'C0') * d2m2R
      )
      d2copContr.d2d2 = d2copContr.d2d2 + (-copContri^-2) * (dDM.d2 %*% t(dDM.d2)) + copContri^-1 * d2DM.d2d2


      # dtheta2.dphi
      d2DM.d2dphi = ( d2copFun.d2dphi(m1R, m2R, phi, type = 'C0') * dm2R
                      + d2copFun.d2dphi(m1L, m2L, phi, type = 'C0') * dm2L
                      - d2copFun.d2dphi(m1R, m2L, phi, type = 'C0') * dm2L
                      - d2copFun.d2dphi(m1L, m2R, phi, type = 'C0') * dm2R
      ) * phi
      d2DM.d2dphi = matrix(d2DM.d2dphi, nrow = length(d2DM.d2dphi))
      d2copContr.d2dphi = d2copContr.d2dphi + (-copContri^-2) * (dDM.d2 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d2dphi


      # dphi.dphi
      d2DM.dphidphi = ( ( d2copFun.dphidphi(m1R, m2R, phi, type = 'C0') + d2copFun.dphidphi(m1L, m2L, phi, type = 'C0')
                          - d2copFun.dphidphi(m1R, m2L, phi, type = 'C0') - d2copFun.dphidphi(m1L, m2R, phi, type = 'C0') ) * phi^2
                        + ( dcopFun.dphi(m1R, m2R, phi, type = 'C0') + dcopFun.dphi(m1L, m2L, phi, type = 'C0')
                            - dcopFun.dphi(m1R, m2L, phi, type = 'C0') - dcopFun.dphi(m1L, m2R, phi, type = 'C0') ) * phi )
      d2copContr.dphidphi = d2copContr.dphidphi + (-copContri^-2) * (dDM.dphi %*% t(dDM.dphi)) + copContri^-1 * d2DM.dphidphi


      rm(m1L); rm(m2L); rm(m1R); rm(m2R); rm(dm1L); rm(dm1R); rm(dm2L); rm(dm2R) # cleaning to avoid invisible bugs
      rm(d2m1L); rm(d2m1R); rm(d2m2L); rm(d2m2R)
      rm(d2DM.d1d1); rm(d2DM.d1d2); rm(d2DM.d1dphi); rm(d2DM.d2d2); rm(d2DM.d2dphi); rm(d2DM.dphidphi)

    } else if(proc1$singleComp$ind.CM$obsT[inds] == 2 & proc2$singleComp$ind.CM$obsT[inds] == 1) {
      m1L = 1 - Ps1[1,nstates]
      m2L = 1 - (Ps2 %*% proc2$singleComp$P.hist[,,jj.marg2])[1,nstates]

      m1R = 1 - (Ps1 %*% proc1$singleComp$P.hist[,,jj.marg1])[1,nstates]

      # LIKELIHOOD
      copContri = copFun(m1L, m2L, phi, type = 'C0') - copFun(m1R, m2L, phi, type = 'C0') # cross contribs and one right time
      copContr = copContr + log(copContri)
      copContrvec = c(copContrvec, copContri)

      # GRADIENT (we will have three vectors here: dtheta1, dtheta2, dphi)
      # dtheta1
      dm1L = -dPs1[1,nstates,]

      dm1R = array(0, dim(proc1$singleComp$dP.hist)[3])
      for(lpar in 1:l.params1) dm1R[lpar] = -(dPs1[,,lpar] %*% proc1$singleComp$P.hist[,,jj.marg1] + Ps1 %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1])[1,nstates]

      dDM.d1 = ( dcopFun.d1(m1L, m2L, phi, type = 'C0') * dm1L - dcopFun.d1(m1R, m2L, phi, type = 'C0') * dm1R )
      dcopContr.d1 = dcopContr.d1 + copContri^-1 * dDM.d1

      # dtheta2
      dm2L = array(0, dim(proc2$singleComp$dP.hist)[3])
      for(lpar in 1:l.params2) dm2L[lpar] = -(dPs2[,,lpar] %*% proc2$singleComp$P.hist[,,jj.marg2] + Ps2 %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2])[1,nstates]

      dDM.d2 = ( dcopFun.d2(m1L, m2L, phi, type = 'C0') * dm2L - dcopFun.d2(m1R, m2L, phi, type = 'C0') * dm2L )
      dcopContr.d2 = dcopContr.d2 + copContri^-1 * dDM.d2

      #dphi
      dDM.dphi = ( dcopFun.dphi(m1L, m2L, phi, type = 'C0') - dcopFun.dphi(m1R, m2L, phi, type = 'C0') ) * phi
      dcopContr.dphi = dcopContr.dphi + copContri^-1 * dDM.dphi


      # HESSIAN

      # (process 1)
      d2m1L = -d2Ps1[1,nstates,,]

      d2m1R = array(0, dim(proc1$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params1){
        for(l2par in 1:l.params1){
          if(l2par >= lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,lpar,l2par,jj.marg1] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,l2par,lpar,jj.marg1]  # lower triangle: take elements from upper triangle and flip them
          d2m1R[lpar, l2par] = -(d2Ps1[,,lpar,l2par] %*% proc1$singleComp$P.hist[,,jj.marg1] + dPs1[,,l2par] %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1]
                                 + dPs1[,,lpar] %*% proc1$singleComp$dP.hist[,,l2par,jj.marg1] + Ps1 %*% d2Pcmpnt)[1,nstates]
        }
      }

      # (process 2)
      d2m2L = array(0, dim(proc2$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params2){
        for(l2par in 1:l.params2){
          if(l2par >= lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,lpar,l2par,jj.marg2] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,l2par,lpar,jj.marg2]  # lower triangle: take elements from upper triangle and flip them
          d2m2L[lpar, l2par] = -(d2Ps2[,,lpar,l2par] %*% proc2$singleComp$P.hist[,,jj.marg2] + dPs2[,,l2par] %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2]
                                 + dPs2[,,lpar] %*% proc2$singleComp$dP.hist[,,l2par,jj.marg2] + Ps2 %*% d2Pcmpnt)[1,nstates]
        }
      }

      # dtheta1.dtheta1
      d2DM.d1d1 = ( d2copFun.d1d1(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm1L)) + dcopFun.d1(m1L, m2L, phi, type = 'C0') * d2m1L
                    - d2copFun.d1d1(m1R, m2L, phi, type = 'C0') * (dm1R %*% t(dm1R)) - dcopFun.d1(m1R, m2L, phi, type = 'C0') * d2m1R
      )
      d2copContr.d1d1 = d2copContr.d1d1 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d1))  + copContri^-1 * d2DM.d1d1

      # dtheta1.dtheta2
      d2DM.d1d2 = ( d2copFun.d1d2(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm2L))
                    - d2copFun.d1d2(m1R, m2L, phi, type = 'C0') * (dm1R %*% t(dm2L))
      )
      d2copContr.d1d2 = d2copContr.d1d2 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d2)) + copContri^-1 * d2DM.d1d2

      # dtheta1.dphi
      d2DM.d1dphi = ( d2copFun.d1dphi(m1L, m2L, phi, type = 'C0') * dm1L
                      - d2copFun.d1dphi(m1R, m2L, phi, type = 'C0') * dm1R
      ) * phi
      d2DM.d1dphi = matrix(d2DM.d1dphi, nrow = length(d2DM.d1dphi))
      d2copContr.d1dphi = d2copContr.d1dphi + (-copContri^-2) * (dDM.d1 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d1dphi

      # dtheta2.dtheta2
      d2DM.d2d2 = ( d2copFun.d2d2(m1L, m2L, phi, type = 'C0') * (dm2L %*% t(dm2L)) + dcopFun.d2(m1L, m2L, phi, type = 'C0') * d2m2L
                    - d2copFun.d2d2(m1R, m2L, phi, type = 'C0') * (dm2L %*% t(dm2L)) - dcopFun.d2(m1R, m2L, phi, type = 'C0') * d2m2L
      )
      d2copContr.d2d2 = d2copContr.d2d2 + (-copContri^-2) * (dDM.d2 %*% t(dDM.d2)) + copContri^-1 * d2DM.d2d2

      # dtheta2.dphi
      d2DM.d2dphi = ( d2copFun.d2dphi(m1L, m2L, phi, type = 'C0') * dm2L
                      - d2copFun.d2dphi(m1R, m2L, phi, type = 'C0') * dm2L
      ) * phi
      d2DM.d2dphi = matrix(d2DM.d2dphi, nrow = length(d2DM.d2dphi))
      d2copContr.d2dphi = d2copContr.d2dphi + (-copContri^-2) * (dDM.d2 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d2dphi

      # dphi.dphi
      d2DM.dphidphi = ( ( d2copFun.dphidphi(m1L, m2L, phi, type = 'C0') - d2copFun.dphidphi(m1R, m2L, phi, type = 'C0') ) * phi^2
                        + ( dcopFun.dphi(m1L, m2L, phi, type = 'C0') - dcopFun.dphi(m1R, m2L, phi, type = 'C0') ) * phi )
      d2copContr.dphidphi = d2copContr.dphidphi + (-copContri^-2) * (dDM.dphi %*% t(dDM.dphi)) + copContri^-1 * d2DM.dphidphi


      rm(m1L); rm(m2L); rm(m1R); rm(dm1L); rm(dm1R); rm(dm2L) # cleaning to avoid invisible bugs
      rm(d2m1L); rm(d2m1R); rm(d2m2L)
      rm(d2DM.d1d1); rm(d2DM.d1d2); rm(d2DM.d1dphi); rm(d2DM.d2d2); rm(d2DM.d2dphi); rm(d2DM.dphidphi)


    } else if(proc1$singleComp$ind.CM$obsT[inds] == 1 & proc2$singleComp$ind.CM$obsT[inds] == 2) {
      m1L = 1 - (Ps1 %*% proc1$singleComp$P.hist[,,jj.marg1])[1,nstates]
      m2L = 1 - Ps2[1,nstates]

      m2R = 1 - (Ps2 %*% proc2$singleComp$P.hist[,,jj.marg2])[1,nstates]

      # LIKELIHOOD
      copContri = copFun(m1L, m2L, phi, type = 'C0') - copFun(m1L, m2R, phi, type = 'C0') # cross contribs and one right time
      copContr = copContr + log(copContri)
      copContrvec = c(copContrvec, copContri)

      # GRADIENT (we will have three vectors here: dtheta1, dtheta2, dphi)
      # dtheta1
      dm1L = array(0, dim(proc1$singleComp$dP.hist)[3])
      for(lpar in 1:l.params1) dm1L[lpar] = -(dPs1[,,lpar] %*% proc1$singleComp$P.hist[,,jj.marg1] + Ps1 %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1])[1,nstates]

      dDM.d1 = ( dcopFun.d1(m1L, m2L, phi, type = 'C0') * dm1L - dcopFun.d1(m1L, m2R, phi, type = 'C0') * dm1L )
      dcopContr.d1 = dcopContr.d1 + copContri^-1 * dDM.d1

      # dtheta2
      dm2L = -dPs2[1,nstates,]

      dm2R = array(0, dim(proc2$singleComp$dP.hist)[3])
      for(lpar in 1:l.params2) dm2R[lpar] = -(dPs2[,,lpar] %*% proc2$singleComp$P.hist[,,jj.marg2] + Ps2 %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2])[1,nstates]

      dDM.d2 = ( dcopFun.d2(m1L, m2L, phi, type = 'C0') * dm2L - dcopFun.d2(m1L, m2R, phi, type = 'C0') * dm2R )
      dcopContr.d2 = dcopContr.d2 + copContri^-1 * dDM.d2

      #dphi
      dDM.dphi = ( dcopFun.dphi(m1L, m2L, phi, type = 'C0') - dcopFun.dphi(m1L, m2R, phi, type = 'C0') ) * phi
      dcopContr.dphi = dcopContr.dphi + copContri^-1 * dDM.dphi

      # HESSIAN

      # (process 1)
      d2m1L = array(0, dim(proc1$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params1){
        for(l2par in 1:l.params1){
          if(l2par >= lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,lpar,l2par,jj.marg1] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,l2par,lpar,jj.marg1]  # lower triangle: take elements from upper triangle and flip them
          d2m1L[lpar, l2par] = -(d2Ps1[,,lpar,l2par] %*% proc1$singleComp$P.hist[,,jj.marg1] + dPs1[,,l2par] %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1]
                                 + dPs1[,,lpar] %*% proc1$singleComp$dP.hist[,,l2par,jj.marg1] + Ps1 %*% d2Pcmpnt)[1,nstates]
        }
      }

      # (process 2)
      d2m2L = -d2Ps2[1,nstates,,]

      d2m2R = array(0, dim(proc2$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params2){
        for(l2par in 1:l.params2){
          if(l2par >= lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,lpar,l2par,jj.marg2] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,l2par,lpar,jj.marg2]  # lower triangle: take elements from upper triangle and flip them
          d2m2R[lpar, l2par] = -(d2Ps2[,,lpar,l2par] %*% proc2$singleComp$P.hist[,,jj.marg2] + dPs2[,,l2par] %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2]
                                 + dPs2[,,lpar] %*% proc2$singleComp$dP.hist[,,l2par,jj.marg2] + Ps2 %*% d2Pcmpnt)[1,nstates]
        }
      }


      # dtheta1.dtheta1
      d2DM.d1d1 = ( d2copFun.d1d1(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm1L)) + dcopFun.d1(m1L, m2L, phi, type = 'C0') * d2m1L
                    - d2copFun.d1d1(m1L, m2R, phi, type = 'C0') * (dm1L %*% t(dm1L)) - dcopFun.d1(m1L, m2R, phi, type = 'C0') * d2m1L
      )
      d2copContr.d1d1 = d2copContr.d1d1 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d1))  + copContri^-1 * d2DM.d1d1

      # dtheta1.dtheta2
      d2DM.d1d2 = ( d2copFun.d1d2(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm2L))
                    - d2copFun.d1d2(m1L, m2R, phi, type = 'C0') * (dm1L %*% t(dm2R))
      )
      d2copContr.d1d2 = d2copContr.d1d2 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d2)) + copContri^-1 * d2DM.d1d2

      # dtheta1.dphi
      d2DM.d1dphi = ( d2copFun.d1dphi(m1L, m2L, phi, type = 'C0') * dm1L
                      - d2copFun.d1dphi(m1L, m2R, phi, type = 'C0') * dm1L
      ) * phi
      d2DM.d1dphi = matrix(d2DM.d1dphi, nrow = length(d2DM.d1dphi))
      d2copContr.d1dphi = d2copContr.d1dphi + (-copContri^-2) * (dDM.d1 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d1dphi

      # dtheta2.dtheta2
      d2DM.d2d2 = ( d2copFun.d2d2(m1L, m2L, phi, type = 'C0') * (dm2L %*% t(dm2L)) + dcopFun.d2(m1L, m2L, phi, type = 'C0') * d2m2L
                    - d2copFun.d2d2(m1L, m2R, phi, type = 'C0') * (dm2R %*% t(dm2R)) - dcopFun.d2(m1L, m2R, phi, type = 'C0') * d2m2R
      )
      d2copContr.d2d2 = d2copContr.d2d2 + (-copContri^-2) * (dDM.d2 %*% t(dDM.d2)) + copContri^-1 * d2DM.d2d2

      # dtheta2.dphi
      d2DM.d2dphi = ( d2copFun.d2dphi(m1L, m2L, phi, type = 'C0') * dm2L
                      - d2copFun.d2dphi(m1L, m2R, phi, type = 'C0') * dm2R
      ) * phi
      d2DM.d2dphi = matrix(d2DM.d2dphi, nrow = length(d2DM.d2dphi))
      d2copContr.d2dphi = d2copContr.d2dphi + (-copContri^-2) * (dDM.d2 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d2dphi

      # dphi.dphi
      d2DM.dphidphi = ( ( d2copFun.dphidphi(m1L, m2L, phi, type = 'C0') - d2copFun.dphidphi(m1L, m2R, phi, type = 'C0') ) * phi^2
                        + ( dcopFun.dphi(m1L, m2L, phi, type = 'C0') - dcopFun.dphi(m1L, m2R, phi, type = 'C0') ) * phi )
      d2copContr.dphidphi = d2copContr.dphidphi + (-copContri^-2) * (dDM.dphi %*% t(dDM.dphi)) + copContri^-1 * d2DM.dphidphi


      rm(m1L); rm(m2L); rm(m2R); rm(dm1L); rm(dm2L); rm(dm2R) # cleaning to avoid invisible bugs
      rm(d2m1L); rm(d2m2L); rm(d2m2R)
      rm(d2DM.d1d1); rm(d2DM.d1d2); rm(d2DM.d1dphi); rm(d2DM.d2d2); rm(d2DM.d2dphi); rm(d2DM.dphidphi)


    } else {
      m1L = 1 - (Ps1 %*% proc1$singleComp$P.hist[,,jj.marg1])[1,nstates]
      m2L = 1 - (Ps2 %*% proc2$singleComp$P.hist[,,jj.marg2])[1,nstates]

      # LIKELIHOOD
      copContri = copFun(m1L, m2L, phi, type = 'C0') # Ps in the two left times
      copContr = copContr + log(copContri)
      copContrvec = c(copContrvec, copContri)

      # GRADIENT (we will have three vectors here: dtheta1, dtheta2, dphi)
      # dtheta1
      dm1L = array(0, dim(proc1$singleComp$dP.hist)[3])
      for(lpar in 1:l.params1) dm1L[lpar] = -(dPs1[,,lpar] %*% proc1$singleComp$P.hist[,,jj.marg1] + Ps1 %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1])[1,nstates]

      dDM.d1 = ( dcopFun.d1(m1L, m2L, phi, type = 'C0') * dm1L )
      dcopContr.d1 = dcopContr.d1 + copContri^-1 * dDM.d1

      # dtheta2
      dm2L = array(0, dim(proc2$singleComp$dP.hist)[3])
      for(lpar in 1:l.params2) dm2L[lpar] = -(dPs2[,,lpar] %*% proc2$singleComp$P.hist[,,jj.marg2] + Ps2 %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2])[1,nstates]

      dDM.d2 = ( dcopFun.d2(m1L, m2L, phi, type = 'C0') * dm2L )
      dcopContr.d2 = dcopContr.d2 + copContri^-1 * dDM.d2

      #dphi
      dDM.dphi = ( dcopFun.dphi(m1L, m2L, phi, type = 'C0') )  * phi
      dcopContr.dphi = dcopContr.dphi + copContri^-1 * dDM.dphi


      # HESSIAN

      # (process 1)
      d2m1L = array(0, dim(proc1$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params1){
        for(l2par in 1:l.params1){
          if(l2par >= lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,lpar,l2par,jj.marg1] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc1$singleComp$d2P.hist[,,l2par,lpar,jj.marg1]  # lower triangle: take elements from upper triangle and flip them
          d2m1L[lpar, l2par] = -(d2Ps1[,,lpar,l2par] %*% proc1$singleComp$P.hist[,,jj.marg1] + dPs1[,,l2par] %*% proc1$singleComp$dP.hist[,,lpar,jj.marg1]
                                 + dPs1[,,lpar] %*% proc1$singleComp$dP.hist[,,l2par,jj.marg1] + Ps1 %*% d2Pcmpnt)[1,nstates]
        }
      }

      # (process 2)
      d2m2L = array(0, dim(proc2$singleComp$d2P.hist)[c(3,4)])
      for(lpar in 1:l.params2){
        for(l2par in 1:l.params2){
          if(l2par >= lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,lpar,l2par,jj.marg2] # upper triangle: OK, we have these in d2P.hist
          if(l2par < lpar) d2Pcmpnt = proc2$singleComp$d2P.hist[,,l2par,lpar,jj.marg2]  # lower triangle: take elements from upper triangle and flip them
          d2m2L[lpar, l2par] = -(d2Ps2[,,lpar,l2par] %*% proc2$singleComp$P.hist[,,jj.marg2] + dPs2[,,l2par] %*% proc2$singleComp$dP.hist[,,lpar,jj.marg2]
                                 + dPs2[,,lpar] %*% proc2$singleComp$dP.hist[,,l2par,jj.marg2] + Ps2 %*% d2Pcmpnt)[1,nstates]
        }
      }


      # dtheta1.dtheta1
      d2DM.d1d1 = ( d2copFun.d1d1(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm1L)) + dcopFun.d1(m1L, m2L, phi, type = 'C0') * d2m1L )
      d2copContr.d1d1 = d2copContr.d1d1 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d1))  + copContri^-1 * d2DM.d1d1

      # dtheta1.dtheta2
      d2DM.d1d2 = ( d2copFun.d1d2(m1L, m2L, phi, type = 'C0') * (dm1L %*% t(dm2L)) )
      d2copContr.d1d2 = d2copContr.d1d2 + (-copContri^-2) * (dDM.d1 %*% t(dDM.d2)) + copContri^-1 * d2DM.d1d2

      # dtheta1.dphi
      d2DM.d1dphi = ( d2copFun.d1dphi(m1L, m2L, phi, type = 'C0') * dm1L ) * phi
      d2DM.d1dphi = matrix(d2DM.d1dphi, nrow = length(d2DM.d1dphi))
      d2copContr.d1dphi = d2copContr.d1dphi + (-copContri^-2) * (dDM.d1 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d1dphi

      # dtheta2.dtheta2
      d2DM.d2d2 = ( d2copFun.d2d2(m1L, m2L, phi, type = 'C0') * (dm2L %*% t(dm2L)) + dcopFun.d2(m1L, m2L, phi, type = 'C0') * d2m2L )
      d2copContr.d2d2 = d2copContr.d2d2 + (-copContri^-2) * (dDM.d2 %*% t(dDM.d2)) + copContri^-1 * d2DM.d2d2

      # dtheta2.dphi
      d2DM.d2dphi = ( d2copFun.d2dphi(m1L, m2L, phi, type = 'C0') * dm2L ) * phi
      d2DM.d2dphi = matrix(d2DM.d2dphi, nrow = length(d2DM.d2dphi))
      d2copContr.d2dphi = d2copContr.d2dphi + (-copContri^-2) * (dDM.d2 %*% t(dDM.dphi)) + copContri^-1 * d2DM.d2dphi

      # dphi.dphi
      d2DM.dphidphi = ( ( d2copFun.dphidphi(m1L, m2L, phi, type = 'C0') ) * phi^2
                        + ( dcopFun.dphi(m1L, m2L, phi, type = 'C0') ) * phi )
      d2copContr.dphidphi = d2copContr.dphidphi + (-copContri^-2) * (dDM.dphi %*% t(dDM.dphi)) + copContri^-1 * d2DM.dphidphi


      rm(m1L); rm(m2L); rm(dm1L); rm(dm2L); # cleaning to avoid invisible bugs
      rm(d2DM.d1d1); rm(d2DM.d1d2); rm(d2DM.d1dphi); rm(d2DM.d2d2); rm(d2DM.d2dphi); rm(d2DM.dphidphi)

    }


  }



  fullDMlik = as.numeric(proc1$singleComp$value + proc2$singleComp$value + (-copContr))
  lDM = as.numeric(proc1$singleComp$l + proc2$singleComp$l + (-copContr))

  # full gradient DM contribution
  fullDMgrad = c(c(proc1$singleComp$gradient) - dcopContr.d1,
                 c(proc2$singleComp$gradient) - dcopContr.d2,
                 -dcopContr.dphi)

  # full Hessian contribution
  DMHess = rbind(cbind(d2copContr.d1d1, d2copContr.d1d2, d2copContr.d1dphi),
                 cbind(t(d2copContr.d1d2), d2copContr.d2d2, d2copContr.d2dphi),
                 cbind(t(d2copContr.d1dphi), t(d2copContr.d2dphi), d2copContr.dphidphi))

  th1.dim = length(params1)
  th2.dim = length(params2)
  phi.dim = length(phi)

  fullDMHess = cbind(rbind(proc1$singleComp$hessian,
                           matrix(0, ncol = th1.dim, nrow = th2.dim+phi.dim)),
                     matrix(0, ncol = th2.dim+phi.dim, nrow = th1.dim+th2.dim+phi.dim)) + cbind(matrix(0, ncol = th1.dim, nrow = th1.dim+th2.dim+phi.dim),
                                                                                                rbind(matrix(0, ncol = th2.dim, nrow = th1.dim),
                                                                                                      proc2$singleComp$hessian, matrix(0, ncol = th2.dim, nrow = phi.dim)),
                                                                                                matrix(0, ncol = phi.dim, nrow = th1.dim+th2.dim+phi.dim)) + (-DMHess)



  # Obtain penalty terms
  S.h  <- cbind(rbind(proc1$singleComp$S.h,
                      matrix(0, ncol = th1.dim, nrow = th2.dim + phi.dim)),
                rbind(matrix(0, ncol = th2.dim, nrow = th1.dim),
                      proc2$singleComp$S.h,
                      matrix(0, ncol = th2.dim, nrow = phi.dim)),
                matrix(0, ncol = phi.dim, nrow = th1.dim + th2.dim + phi.dim)) # hess
  S.h1 <- 0.5*crossprod(params, S.h) %*% params                                # lik
  S.h2 <- S.h %*% params                                                       # grad



  print(paste(Sys.time(), ' // Lik = ', c(fullDMlik), ' // Max abs grad = ',
              max(abs(fullDMgrad)), ' // min eigen Hess = ',
              min(eigen(fullDMHess)$values)))


  list(value = fullDMlik,
       gradient = fullDMgrad,
       hessian = fullDMHess,
       S.h = S.h,
       S.h1 = S.h1,
       S.h2 = S.h2,
       l = lDM,
       proc1 = proc1,
       proc2 = proc2)


}
