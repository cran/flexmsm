

#' Likelihood, gradient and Hessian for univariate transition intensity based models for the base dependence model.
#'
#' @param params Parameters vector.
#' @param data Dataset in proper format.
#' @param full.X Full design matrix.
#' @param MM List of necessary setup quantities.
#' @param pen.matr.S.lambda Penalty matrix multiplied by smoothing parameter lambda.
#' @param aggregated.provided Whether aggregated form was provided (may become obsolete in the future if we see original dataset as special case of aggregated where \code{nrep = 1}).
#' @param do.gradient Whether or not to compute the gradient.
#' @param do.hessian Whether or not to compute the Hessian.
#' @param pmethod Method to be used for computation of transition probability matrix. See help of \code{msm()} for further details.
#' @param death Whether the last state is an absorbing state.
#' @param Qmatr.diagnostics.list List of maximum absolute values of the Q matrices computed during model fitting.
#' @param verbose Whether to print out the progress being made in computing the likelihood, gradient and Hessian.
#' @param parallel Whether or not to use parallel computing (only for Windows users for now).
#' @param no_cores Number of cores used if parallel computing chosen. The default is 2. If \code{NULL}, all available cores are used.
#' @param CM.comp If \code{TRUE} compute the Copula contribution from each individual for the base dependence model. The default is \code{TRUE}.
#' @param P.save.all If \code{TRUE} save the P matrices computed at each observation. The default is \code{FALSE}.
#'
#' @return Penalized likelihood, gradient and Hessian associated with model at given parameters, for use by trust region algorithm.
#' @export
#'
#'
LikGradHess.CM = function(params, data = NULL, full.X = NULL, MM, pen.matr.S.lambda, aggregated.provided = FALSE,
                          do.gradient = TRUE, do.hessian = TRUE, pmethod = 'analytic', death,
                          Qmatr.diagnostics.list = NULL,
                          verbose = FALSE, parallel = FALSE, no_cores = 2,
                          CM.comp = TRUE, P.save.all = FALSE){



  # Initialize
  l.par = 0                                             # log-likelihood
  G = rep(0, MM$l.params)                               # gradient
  H = matrix(0, ncol = MM$l.params, nrow = MM$l.params) # hessian
  apprH = matrix(0, ncol = MM$l.params, nrow = MM$l.params) # approximate hessian

  if(is.null(MM$cens.state)) MM$cens.state = -99

  P.hist = dP.hist = d2P.hist = NULL
  if(P.save.all) {
    i.save = 1 # for the saving of the P matrix, if required
    P.hist = array(0, dim = c(MM$nstate, MM$nstates, nrow(data)-length(unique(data$"(id)"))))   # saving the P matrices if required
    dP.hist = array(0, dim = c(MM$nstate, MM$nstates, MM$l.params, nrow(data)-length(unique(data$"(id)"))))  # saving the dP matrices if required
    d2P.hist = array(0, dim = c(MM$nstate, MM$nstates, MM$l.params, MM$l.params, nrow(data)-length(unique(data$"(id)")))) # saving the d2P matrices if required
  }

  # this is for the contributions to the dependence model
  P.CM.contr = list(list(diag(nrow = MM$nstates)))

  # indicator joining point
  ind.CM = data.frame(v1 = data$"(id)"[1], obsT = 1, ni = sum(data$"(id)" == data$"(id)"[1])-1)
  names(ind.CM)[1] = "(id)"

  if(parallel){
    if(is.null(no_cores)) no_cores <- parallel::detectCores(logical = TRUE)  # returns the number of available hardware threads, and if it is FALSE, returns the number of physical cores
    if(!is.null(no_cores) & no_cores > parallel::detectCores(logical = TRUE)) stop('You can use up to ', parallel::detectCores(logical = TRUE), ' cores on this device. Please change the value provided for no_cores.')

    # Break the data in no_cores sections
    peeps = unique(data$`(id)`)
    part.lik = floor(length(peeps)/no_cores)
    part.lik.last = length(peeps)  - part.lik*(no_cores-1)

    sta = 1; sto = sta + part.lik - 1
    data.list = list(list(data[data$`(id)` %in% peeps[sta:sto], ]))
    data.list[[1]][[2]] = full.X[which(data$`(id)` %in% peeps[sta:sto]), ]

    for(i in 2:no_cores){
      sta = sto + 1
      sto = sta + ifelse(i < no_cores, part.lik, part.lik.last) - 1
      data.list[[i]] = list(data[data$`(id)` %in% peeps[sta:sto], ])
      data.list[[i]][[2]] = full.X[which(data$`(id)` %in% peeps[sta:sto]), ]
    }

    parallel.LikGradHess = function(data.parallel) {
      LikGradHess.CM(params, data = data.parallel[[1]], full.X = data.parallel[[2]], MM = MM, pen.matr.S.lambda = matrix(0, nrow = nrow(pen.matr.S.lambda), ncol = ncol(pen.matr.S.lambda)), # *****
                          aggregated.provided = FALSE, do.gradient = do.gradient, do.hessian = do.hessian,
                          pmethod = pmethod, death = death,
                          Qmatr.diagnostics.list = Qmatr.diagnostics.list,
                          parallel = FALSE, P.save.all = P.save.all, CM.comp = CM.comp)
    }


    cl <- parallel::makeCluster(no_cores)
    # parallel::clusterEvalQ(cl, {library(flexmsm)})

    parallel::clusterExport(cl, varlist = c('LikGradHess.CM', 'Q.matr.setup.general', 'P.matr.comp', 'dP.matr.comp', 'd2P.matr.comp',
                                            'params', 'data', 'full.X', 'MM', 'pen.matr.S.lambda',
                                            'aggregated.provided', 'do.gradient', 'do.hessian',
                                            'pmethod', 'death',
                                            'Qmatr.diagnostics.list', 'P.save.all', 'CM.comp'), envir = environment())
    parallel.LikGradHess <- parallel::parLapply(cl, data.list, parallel.LikGradHess)
    parallel::stopCluster(cl)



    Q.diag.list.temp = c()
    jj.cores.tmp = 1
    ind.CM = parallel.LikGradHess[[1]]$ind.CM
    P.CM.contr = parallel.LikGradHess[[1]]$P.CM.contr
    for(i in 1:no_cores){ # bring together sub-quantities computed at each core
      l.par = l.par + parallel.LikGradHess[[i]]$value
      G = G + parallel.LikGradHess[[i]]$gradient
      H = H + parallel.LikGradHess[[i]]$hessian
      Q.diag.list.temp = c(Q.diag.list.temp, parallel.LikGradHess[[i]]$Qmatr.diagnostics.list[[length(parallel.LikGradHess[[i]]$Qmatr.diagnostics.list)]])
      if(P.save.all){
        P.hist[,,jj.cores.tmp:(jj.cores.tmp + dim(parallel.LikGradHess[[i]]$P.hist)[3]-1)] = parallel.LikGradHess[[i]]$P.hist
        dP.hist[,,,jj.cores.tmp:(jj.cores.tmp + dim(parallel.LikGradHess[[i]]$dP.hist)[4]-1)] = parallel.LikGradHess[[i]]$dP.hist
        d2P.hist[,,,,jj.cores.tmp:(jj.cores.tmp + dim(parallel.LikGradHess[[i]]$d2P.hist)[5]-1)] = parallel.LikGradHess[[i]]$d2P.hist
        if(i > 1) ind.CM = rbind(ind.CM, parallel.LikGradHess[[i]]$ind.CM)

        jj.cores.tmp = jj.cores.tmp + dim(parallel.LikGradHess[[i]]$P.hist)[3]
      }
      if(CM.comp) if(i > 1) for(l.inn in 1:length(parallel.LikGradHess[[i]]$P.CM.contr)) P.CM.contr[[length(P.CM.contr)+1]] = parallel.LikGradHess[[i]]$P.CM.contr[[l.inn]]
    }
    Qmatr.diagnostics.list[[length(Qmatr.diagnostics.list)+1]] = Q.diag.list.temp
    rm(jj.cores.tmp)

    l.par = -l.par; G = -G; H = -H # because we will be flipping the sign later, when we add the real final penalty



  } else { # **** ORIGINAL (NOT PARALLELISED) CODE ****



    Q.matr.object = Q.matr.setup.general(params[MM$pos.optparams], MM$nstates, full.X, MM$start.pos.par, MM$l.short.formula, MM$whereQ,
                                         firstD = do.gradient, secondD = do.hessian, bound.eta = FALSE, # note: bound.eta was only for debugging purposes, don't touch this
                                         pos.optparams = MM$pos.optparams, pos.optparams2 = MM$pos.optparams2)
    Qmatr = Q.matr.object$Qmatr
    dQmatr = Q.matr.object$dQmatr
    d2Qmatr = Q.matr.object$d2Qmatr
    comp.par.mapping = Q.matr.object$comp.par.mapping


    # [03/01/2021] Introducing a check to verify whether pathological behaviour is occurring within the Q matrix and when this starts
    # to happen ------------------------------------------------------------------------- PART A -----------------------------------
    if(!is.null(Qmatr.diagnostics.list)) Qmatr.i.maxs = c()
    # ------------------------------------------------------------------------------------------------------------------------------

    # # Set first id, Q and dQ (then only fish out new Q matrix and new dQ matrix when this changes) --- only id
    # id.t = data$"(id)"[1]
    # # *************************************************************** NOTE: id.t is never actually used, remove tentatively


    for(i in 1:nrow(data)){ # perhaps change to while(i < nagg) where nagg is the number of rows in the aggregated dataset as done in Jackson lik.c code (on github)

      # print(i)
      # # For checking the speed of this new implementation
      # print(paste(i, '/'  , nrow(data), sep = ''))

      if( i != nrow(data) & data$"(id)"[i] != data$"(id)"[i+1] ){
        P.CM.contr[[length(P.CM.contr)+1]] = list(diag(nrow = MM$nstates)) # create the new CM contribution when the individual changes
        ind.CM[nrow(ind.CM)+1,] = c(data$"(id)"[i+1], 1, sum(data$"(id)" == data$"(id)"[i+1])-1)
      }

      if( (i == nrow(data) || data$"(id)"[i] != data$"(id)"[i+1]) ) next # i.e. skip the last observation for each individual

      # if( data$"(id)"[i] != id.t ) id.t = data$"(id)"[i] # update id.t when it changes

      fromstate = data$"(fromstate)"[i]
      tostate = data$"(tostate)"[i]
      timelag = data$"(timelag)"[i]
      nrep = data$"(nrep)"[i]
      living.exact = data$"(living.exact)"[i]

      Qmatr.i   = Qmatr[,,i]
      dQmatr.i  = dQmatr[,,,i]
      d2Qmatr.i = d2Qmatr[,,,i] #d2Qmatr[,,,,i]

      # [03/01/2021] Introducing a check to verify whether pathological behaviour is occurring within the Q matrix and when this starts
      # to happen ------------------------------------------------------------------------- PART B -----------------------------------

      if(!is.null(Qmatr.diagnostics.list)) Qmatr.i.maxs = c(Qmatr.i.maxs,  max(abs(Qmatr.i)))
      # --------------------------------------------------------------------------------------------------------------------------------

      if(pmethod == 'eigendecomp'){
        # *** EIGENDECOMPISITION ***
        # Compute P matrix and extract eigendecomposition quantities
        # times.eig = c()
        # start.eig = Sys.time()
        P.comp.i = P.matr.comp(Qmatr.i, timelag)

        Pmatr.i = P.comp.i$Pmatr
        decomp.obj.i = P.comp.i$decomp.obj

        # Compute dP matrices (this is an array of dimension nstates x nstates x l.params)
        if(do.gradient) dPmatr.i = dP.matr.comp(nstates = MM$nstates, l.params = MM$l.params, timelag = timelag,
                                                dQmatr.i = dQmatr.i, method = 'eigendecomp', decomp.obj.i = decomp.obj.i)$dPmatr

        # Compute d2P matrices (this is an array of dimension nstates x nstates x l.params x l.params)
        if(do.hessian) d2Pmatr.i = d2P.matr.comp(nstates = MM$nstates, l.params = MM$l.params, timelag, dQmatr.i, d2Qmatr.i, method = 'eigendecomp',
                                                 decomp.obj.i = decomp.obj.i, start.pos.par = MM$start.pos.par, comp.par.mapping = comp.par.mapping)

      }



      if(pmethod == 'analytic'){ # DO WE NEED COMP.PAR.MAPPING HERE TOO?? PROBABLY! REMEMBER TO CHECK THIS!!
        # *** ANALYTIC EXPRESSIONS ***
        start.anal = Sys.time()
        Pmatr.i = P.matr.comp(Qmatr = Qmatr.i, timelag = timelag, method = 'analytic')$Pmatr

        dPmatr.i = dP.matr.comp(nstates = MM$nstates, l.params = MM$l.params, timelag = timelag,
                                dQmatr.i = dQmatr.i, method = 'analytic',
                                Pmatr.i = Pmatr.i, Qmatr.i = Qmatr.i)$dPmatr

        d2Pmatr.i = d2P.matr.comp(nstates = MM$nstates, l.params = MM$l.params, timelag = timelag,
                                  dQmatr.i = dQmatr.i, d2Qmatr.i = d2Qmatr.i, start.pos.par = MM$start.pos.par,
                                  method = 'analytic',
                                  Qmatr.i = Qmatr.i, Pmatr.i = Pmatr.i, dPmatr.i = dPmatr.i)

        end.anal = Sys.time()
      }



      if(pmethod == 'scaling&squaring'){
        # *** SCALING AND SQUARING ***
        start.new = Sys.time()
        Pmatr.comp.i = P.matr.comp(Qmatr = Qmatr.i, timelag = timelag, method = 'scaling&squaring')

        Pmatr.i =  Pmatr.comp.i$Pmatr
        SS.obj.i = Pmatr.comp.i$SS.obj

        dPmatr.comp.i = dP.matr.comp(nstates = MM$nstates, l.params = MM$l.params, timelag = timelag,
                                     dQmatr.i = dQmatr.i, method = 'scaling&squaring', Qmatr.i = Qmatr.i, SS.obj.i = SS.obj.i)
        dPmatr.i.SS = dPmatr.comp.i$dPmatr
        SS.obj.i = dPmatr.comp.i$SS.obj

        d2Pmatr.i.SS = d2P.matr.comp(nstates = MM$nstates, l.params = MM$l.params, timelag = timelag,
                                     dQmatr.i = dQmatr.i, d2Qmatr.i = d2Qmatr.i, start.pos.par = MM$start.pos.par,
                                     method = 'scaling&squaring',
                                     Qmatr.i = Qmatr.i, SS.obj.i = SS.obj.i)
        end.new = Sys.time()

        dPmatr.i = dPmatr.i.SS
        d2Pmatr.i = d2Pmatr.i.SS
      }

      # save the P matrix if required
      if(P.save.all){
        P.hist[,,i.save] = Pmatr.i
        if(do.gradient) dP.hist[,,,i.save] = dPmatr.i
        if(do.hessian) d2P.hist[,,,,i.save] = d2Pmatr.i
        i.save = i.save + 1
      }

      # DEPENDENCE MODEL CONTRIBUTION
      if(CM.comp){
        if(tostate != MM$nstates){
          P.CM.contr[[length(P.CM.contr)]][[1]] = P.CM.contr[[length(P.CM.contr)]][[1]] %*% Pmatr.i
        }

        if(tostate == MM$nstates & fromstate != MM$nstates){
          P.CM.contr[[length(P.CM.contr)]][[2]] = P.CM.contr[[length(P.CM.contr)]][[1]] %*% Pmatr.i
          ind.CM$obsT[nrow(ind.CM)] = 2
        }
      }



      # do.gradient = TRUE # might be useful for debugging purposes
      # do.hessian = TRUE # might be useful for debugging purposes

      if((tostate != MM$nstates & fromstate != MM$cens.state & tostate != MM$cens.state & !living.exact) | !death) {

        # log-likelihood --------------------------
        l.par = l.par + log(Pmatr.i[fromstate, tostate]) * nrep

        # These are just checks used in debug
        if( is.nan(log(Pmatr.i[fromstate, tostate])) ) stop(paste('A NaN was produced at observation', i))
        if( -log(Pmatr.i[fromstate, tostate]) == Inf ) stop(paste('A Inf was produced at observation', i))
        # print(paste('i =', i, log(Pmatr.i[fromstate, tostate]) ))



        if(do.gradient){

          # gradient --------------------------------------------------------------------------------------
          G = G + dPmatr.i[fromstate, tostate, ]/Pmatr.i[fromstate, tostate] * nrep # exploiting array structure

        }


        if(do.hessian){

          # hessian -----------------------------------------------------------------------------------------
          second.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
          prod.first.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
          for(k in 1:MM$l.params){
            for(l in k:MM$l.params){
              prod.first.deriv.Ind.i[k,l] = dPmatr.i[fromstate, tostate, k] * dPmatr.i[fromstate, tostate, l]
              second.deriv.Ind.i[k,l] = d2Pmatr.i[fromstate, tostate, k, l]
            }
          }
          H = H + (second.deriv.Ind.i/Pmatr.i[fromstate, tostate] - prod.first.deriv.Ind.i/Pmatr.i[fromstate, tostate]^2) * nrep


          # Approximate Hessian -----------------------------------------------------------------------------------------
          apprH = apprH + (#second.deriv.Ind.i/Pmatr.i[fromstate, tostate]
            - prod.first.deriv.Ind.i/Pmatr.i[fromstate, tostate]^2) * nrep

        }
      }
      # else if (tostate == MM$nstates & death) { # absorbing state
      #
      #   # log-likelihood -------------------------------------------------------------------------------------------
      #   exact.death.term = 0
      #   for(s in 1:(MM$nstates-1)) exact.death.term = exact.death.term + Pmatr.i[fromstate, s] * Qmatr.i[s, MM$nstates]
      #   l.par = l.par + log(exact.death.term) * nrep
      #
      #   if(do.gradient){
      #
      #     # gradient -------------------------------------------------------------------------------------------------
      #     exact.death.term = 0
      #     exact.death.term.der = 0
      #     for(s in 1:(MM$nstates-1)){
      #       exact.death.term = exact.death.term + Pmatr.i[fromstate, s]*Qmatr.i[s, MM$nstates]
      #       exact.death.term.der = exact.death.term.der + (dPmatr.i[fromstate, s, ]*Qmatr.i[s, MM$nstates] + Pmatr.i[fromstate, s]*dQmatr.i[s, MM$nstates, ])
      #     }
      #
      #     G = G + exact.death.term.der/exact.death.term * nrep
      #
      #   }
      #
      #   if(do.hessian){
      #
      #     # hessian -----------------------------------------------------------------------------------------------------
      #     exact.death.term = 0
      #     single.term.1 = 0
      #     single.term.2 = 0
      #     square.term = 0
      #
      #
      #
      #     for(s in 1:(MM$nstates-1)){
      #
      #       first.deriv.term.1 = rep(0, MM$l.params)
      #       first.deriv.term.2 = matrix(0, MM$l.params)
      #       second.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #       second.deriv.Q.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #       prod.first.deriv.Ind.i.1 = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #       prod.first.deriv.Ind.i.2 = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #
      #       comp.par = 1 # to properly map d2Q given its compact form
      #
      #       for(k in 1:MM$l.params){
      #
      #         first.deriv.term.1[k] = dPmatr.i[fromstate, s, k]*Qmatr.i[s, MM$nstates] + Pmatr.i[fromstate, s]*dQmatr.i[s, MM$nstates, k]
      #
      #         for(l in k:MM$l.params){ # only computation of upper triangle is needed (reconstruct full matrix from this later)
      #
      #           prod.first.deriv.Ind.i.1[k,l] = dPmatr.i[fromstate, s, k]*dQmatr.i[s, MM$nstates, l]
      #           prod.first.deriv.Ind.i.2[k,l] = dQmatr.i[s, MM$nstates, k]*dPmatr.i[fromstate, s, l]
      #
      #           second.deriv.Ind.i[k,l] = d2Pmatr.i[fromstate, s, k, l]
      #
      #           # COMPACT IMPLEMENTATION (remember to uncomment comp.par as well) **************
      #           if( is.null(comp.par.mapping) ){
      #             if( max(which(k - MM$start.pos.par >= 0)) == max(which(l - MM$start.pos.par >= 0)) ) { # ... but remember block-diagonal structure, this checks to which transition intensity the parameter belongs, if not same for r1 and r2 then d2Q is just matrix of zeroes
      #               second.deriv.Q.Ind.i[k,l] = d2Qmatr.i[s, MM$nstates, comp.par]
      #               comp.par = comp.par + 1
      #             } else {
      #               second.deriv.Q.Ind.i[k,l] = 0
      #             }
      #           } else {
      #
      #             if( !is.na(comp.par.mapping[k, l]) ){
      #               second.deriv.Q.Ind.i[k,l] = d2Qmatr.i[s, MM$nstates, comp.par.mapping[k, l]]
      #             } else {
      #               second.deriv.Q.Ind.i[k,l] = 0
      #             }
      #
      #           }
      #           # *********************************************************************************
      #
      #         }
      #       }
      #
      #       exact.death.term = exact.death.term + Pmatr.i[fromstate, s]*Qmatr.i[s, MM$nstates]
      #       single.term.1 = single.term.1 + first.deriv.term.1
      #       square.term = square.term + (second.deriv.Ind.i*Qmatr.i[s, MM$nstates] + prod.first.deriv.Ind.i.1 + prod.first.deriv.Ind.i.2 + Pmatr.i[fromstate, s]*second.deriv.Q.Ind.i)
      #
      #     }
      #     matr.single.term.1 = matrix(rep(single.term.1, MM$l.params), ncol = MM$l.params, byrow = TRUE)
      #     matr.single.term.2 = t(matr.single.term.1)
      #
      #     matr.single.term.1[lower.tri(matr.single.term.1, diag = FALSE)] = 0 # only upper triangle
      #     matr.single.term.2[lower.tri(matr.single.term.2, diag = FALSE)] = 0 # only upper triangle
      #
      #     H = H + (square.term/exact.death.term - matr.single.term.1*matr.single.term.2/exact.death.term^2) * nrep
      #
      #   }
      #
      #
      # } else if (living.exact) { # exactly observed living state
      #
      #
      #   l.par = l.par + ( Qmatr.i[fromstate, fromstate] * timelag + log(Qmatr[fromstate, tostate, i+1]) ) * nrep # last term is 'i+1' because the Q we multiply by is the one found in the arrival time, and since we are assuming left-cont piecewise const we need the Q from the following observation
      #
      #   if(do.gradient) G = G + (dQmatr[fromstate, tostate, , i+1] -  Qmatr[fromstate, tostate, i+1] * dQmatr.i[fromstate, fromstate, ] *  timelag ) / Qmatr[fromstate, tostate, i+1] * nrep
      #
      #   if(do.hessian){
      #
      #     prod.first.deriv.Ind.i.1 = prod.first.deriv.Ind.i.2 = prod.first.deriv.Ind.i.3 = prod.first.deriv.Ind.i.3 = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     prod.first.deriv.Ind.i.4 = prod.first.deriv.Ind.i.5 = prod.first.deriv.Ind.i.6 = prod.first.deriv.Ind.i.7 = prod.first.deriv.Ind.i.3 = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     second.deriv.Q.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     second.deriv.Q.Ind.ip1 = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #
      #     comp.par = 1
      #
      #     for(k in 1:MM$l.params){
      #       for (l in k:MM$l.params){
      #         prod.first.deriv.Ind.i.1[k,l] =  dQmatr[fromstate, tostate, k, i+1] * dQmatr.i[fromstate, fromstate, l] * timelag
      #         prod.first.deriv.Ind.i.2[k,l] =  Qmatr[fromstate, tostate, k, i+1] * dQmatr.i[fromstate, fromstate, k] * dQmatr.i[fromstate, fromstate, l] * timelag^2
      #         prod.first.deriv.Ind.i.3[k,l] =  dQmatr[fromstate, tostate, l, i+1] * dQmatr.i[fromstate, fromstate, k] * timelag
      #
      #         prod.first.deriv.Ind.i.4[k,l] =  dQmatr[fromstate, tostate, k, i+1] * dQmatr[fromstate, tostate, l, i+1]
      #         prod.first.deriv.Ind.i.5[k,l] =  Qmatr[fromstate, tostate, i+1] * dQmatr[fromstate, tostate, k, i+1] * dQmatr.i[fromstate, fromstate, l] * timelag
      #         prod.first.deriv.Ind.i.6[k,l] =  Qmatr[fromstate, tostate, i+1] * dQmatr.i[fromstate, tostate, k] * dQmatr[fromstate, fromstate, l, i+1] * timelag
      #         prod.first.deriv.Ind.i.7[k,l] =  Qmatr[fromstate, tostate, i+1]^2 * dQmatr.i[fromstate, tostate, k] * dQmatr.i[fromstate, tostate, l] * timelag^2
      #
      #         # COMPACT IMPLEMENTATION (remember to uncomment comp.par as well) **************
      #         if( is.null(comp.par.mapping) ){
      #           if( max(which(k - MM$start.pos.par >= 0)) == max(which(l - MM$start.pos.par >= 0)) ) { # ... but remember block-diagonal structure, this checks to which transition intensity the parameter belongs, if not same for r1 and r2 then d2Q is just matrix of zeroes
      #             second.deriv.Q.Ind.i[k,l] = d2Qmatr.i[fromstate, fromstate, comp.par]
      #             second.deriv.Q.Ind.ip1[k,l] = d2Qmatr[fromstate, tostate, comp.par,i+1]
      #
      #             comp.par = comp.par + 1
      #           } else {
      #             second.deriv.Q.Ind.i[k,l] = 0
      #           }
      #         } else {
      #
      #           if( !is.na(comp.par.mapping[k, l]) ){
      #             second.deriv.Q.Ind.i[k,l] =   d2Qmatr.i[fromstate, fromstate, comp.par.mapping[k, l]]
      #             second.deriv.Q.Ind.ip1[k,l] = d2Qmatr[fromstate, tostate, comp.par.mapping[k, l],i+1]
      #
      #           } else {
      #             second.deriv.Q.Ind.i[k,l] = 0
      #             second.deriv.Q.Ind.ip1[k,l] = 0
      #           }
      #
      #         }
      #         # *********************************************************************************
      #       }
      #     }
      #
      #     H = H + (Qmatr[fromstate, tostate, i+1]^-1 * ( prod.first.deriv.Ind.i.1 + prod.first.deriv.Ind.i.2 + second.deriv.Q.Ind.ip1 + prod.first.deriv.Ind.i.3 + Qmatr[fromstate, tostate, i+1] * second.deriv.Q.Ind.i * timelag) + Qmatr[fromstate, tostate, i+1]^-2 * ( prod.first.deriv.Ind.i.4 + prod.first.deriv.Ind.i.5 + prod.first.deriv.Ind.i.6 + prod.first.deriv.Ind.i.7 ) ) * nrep
      #
      #   }
      #
      # } else if (tostate == MM$cens.state) { # censored state
      #
      #   # log-likelihood --------------------------
      #   l.par = l.par + log( sum(Pmatr.i[fromstate, ]) ) * nrep
      #
      #   # gradient --------------------------------------------------------------------------------------
      #   if(do.gradient) G = G + apply(dPmatr.i[fromstate, , ], MARGIN = 2, FUN = sum)/sum(Pmatr.i[fromstate, ]) * nrep # exploiting array structure
      #
      #
      #
      #   if(do.hessian){
      #
      #     # hessian -----------------------------------------------------------------------------------------
      #     second.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     prod.first.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     for(k in 1:MM$l.params){
      #       for(l in k:MM$l.params){
      #         prod.first.deriv.Ind.i[k,l] = sum(dPmatr.i[fromstate, , k]) * sum(dPmatr.i[fromstate, , l])
      #         second.deriv.Ind.i[k,l] = sum(d2Pmatr.i[fromstate, , k, l])
      #       }
      #     }
      #     H = H + (second.deriv.Ind.i/sum(Pmatr.i[fromstate, ]) - prod.first.deriv.Ind.i/sum(Pmatr.i[fromstate, ])^2) * nrep
      #
      #   }
      #
      # } else if (fromstate == MM$cens.state) { # censored state starting point - not sure if this is correct !! keep an eye on it
      #
      #   # log-likelihood --------------------------
      #   l.par = l.par + log( sum(Pmatr.i[, tostate]) ) * nrep
      #
      #   # gradient --------------------------------------------------------------------------------------
      #   if(do.gradient) G = G + apply(dPmatr.i[, tostate, ], MARGIN = 2, FUN = sum)/sum(Pmatr.i[, tostate]) * nrep # exploiting array structure
      #
      #
      #
      #   if(do.hessian){
      #
      #     # hessian -----------------------------------------------------------------------------------------
      #     second.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     prod.first.deriv.Ind.i = matrix(0, ncol = MM$l.params, nrow = MM$l.params)
      #     for(k in 1:MM$l.params){
      #       for(l in k:MM$l.params){
      #         prod.first.deriv.Ind.i[k,l] = sum(dPmatr.i[, tostate, k]) * sum(dPmatr.i[, tostate, l])
      #         second.deriv.Ind.i[k,l] = sum(d2Pmatr.i[, tostate, k, l])
      #       }
      #     }
      #     H = H + (second.deriv.Ind.i/sum(Pmatr.i[, tostate]) - prod.first.deriv.Ind.i/sum(Pmatr.i[, tostate])^2) * nrep
      #
      #   }
      #
      # }
    }



  }



  # Obtain penalty terms
  S.h  <- pen.matr.S.lambda                      # hess
  S.h1 <- 0.5*crossprod(params, S.h)%*%params    # lik
  S.h2 <- S.h%*%params                           # grad

  # Penalize log-likelihood
  lik.pen = -l.par + S.h1

  # Penalize gradient
  G.pen = -G + S.h2

  # Penalize hessian
  # H[lower.tri(H, diag = FALSE)] = H[upper.tri(H, diag = FALSE)] # reconstruct full matrix --- THIS IS WRONG !!!
  if(!parallel){ # no need as we already return the full matrix with parallel
    H = H + t(H)         # reconstruct full matrix
    diag(H) = diag(H)/2  # reconstruct full matrix

    apprH = apprH + t(apprH)         # reconstruct full matrix
    diag(apprH) = diag(apprH)/2  # reconstruct full matrix
  }
  H.pen = -H + S.h
  apprH.pen = -apprH + S.h


  # Q DIAGNOTICS - may be useful
  if(!is.null(Qmatr.diagnostics.list) & !parallel){
    Qmatr.diagnostics.list[[length(Qmatr.diagnostics.list)+1]] = Qmatr.i.maxs
  }

  # --------------------------------------------------------------------------------------------------------------------------------


  if(verbose){
    cat("One likelihood, gradient and Hessian call fully completed!", Sys.time())
    cat('lik.pen =', lik.pen, "// max abs grad =", round(max(abs(G.pen)), 5), '// min eigen =', round(min(eigen(H.pen)$values, 5)))
  }


  list(value = lik.pen, gradient = G.pen, hessian = H.pen, # penalized log-lik, gradient and Hessian
       apprHessian = apprH,                                # Hessian obtained with outer product approximation
       S.h = S.h, S.h1 = S.h1, S.h2 = S.h2,                # penalty matrices
       l = -l.par,                                         # unpenalized negative log-likelihood
       Qmatr.diagnostics.list = Qmatr.diagnostics.list,    # containing maximum abs Q values, for diagnostics
       P.hist = P.hist,                                    # saving the P matrix for each observation if required
       dP.hist = dP.hist,                                  # saving the dP matrix for each observation if required
       d2P.hist = d2P.hist,                                # saving the d2P matrix for each observation if required
       P.CM.contr = P.CM.contr,                            # saving only the contributions needed for the CM dependence structure
       ind.CM = ind.CM                                     # saving indicator to be used for the computation of the P, dP and d2P terms of the CM
  )

}
