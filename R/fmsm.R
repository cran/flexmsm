#' Flexible transition intensity based models for univariate multistate processes
#'
#' @description
#' Fits a flexible multistate survival model. Any type of process is supported,
#' including both forward and backward transitions, and must be specified by providing a list
#' of equations, one for each transition intensity allowed. Any type of observation scheme is
#' allowed: the process can be observed in continuous time, intermittently at fixed times, there
#' can be an absorbing state as well as censored states.
#'  Virtually any type of covariate
#' effects are supported and can be specified by means of splines, with the same syntax used
#' to specify Generalised Additive Models (GAMs) in \code{R}.
#'
#' @param formula Model specification for the transition intensities.
#' @param data Dataset.
#' @param id Name of the variable in the dataset representing the unique code associated with each patient.
#' @param state Name of the variable in the dataset representing the state occupied by the patient at the given time.
#' @param death \code{TRUE} if the last state is an absorbing state, \code{FALSE} otherwise.
#' @param pmethod Which method should be used for the computation of the transition probability matrix. Available options are
#' \itemize{
#' \item \code{'eigendecomp'} (default): this method is based on the eigendecomposition of the transition intensity matrix (from Kalbfleisch & Lawless 1985);
#' \item \code{'analytic'}: uses analytic expressions of the transition probabilities, obtained by solving the Kolmogorov forward differential equation, only implemented for IDMs for now;
#' \item \code{'scaling&squaring'}: this is the scaling and squaring method implemented as proposed in Fung (2004).This is inefficient, so its use is not recommended. Can be used to investigate convergence errors.
#' }
#' @param aggregate Whether or not data should be aggregated (this slightly improves efficiency as redundancies in the data are eliminated). The default is \code{TRUE}.
#' @param params.0 Starting values for the model parameters. Defaults to \code{NULL}, i.e. they are computed internally.
#' @param sp.0 Starting values for the smoothing parameters. Defaults to \code{NULL}, i.e. they are computed internally.
#' @param constraint A list containing the constraints to be applied to the model parameters. For example, assuming a process with three transitions, \code{constraint = list(x1 = c(1,1,1), x2 = c(1,2,2))}  means that the effect of covariate \code{x1} is constrained to be equal across the tree transitions and that the effects of \code{x2} on the second and third transitions are constrained to be equal, but the effect on the first transition is left unconstrained.
#' @param sp.method Method to be used for smoothing parameter estimation. The default is \code{magic}, the automatic multiple smoothing parameter selection algorithm. Alternatively, \code{efs} can be used for the Fellner-Schall method. To suppress the smoothing parameter estimation set this to \code{NULL}.
#' @param iterlimsp Maximum allowed iterations for smoothing parameter estimation.
#' @param Q.diagnostics If \code{TRUE}, diagnostics information on the Q matrix are saved. The default \code{TRUE}.
#' @param fit If \code{FALSE}, fitting is not carried. May be useful to extract model setup quantities.
#' @param iterlim Maximum allowed iterations for trust region algorithm.
#' @param tolsp Convergence criterion used in \code{magic} based smoothing parameter estimation.
#' @param tolsp.EFS Convergence criterion used in \code{efs} based smoothing parameter estimation.
#' @param parallel If \code{TRUE} parallel computing is used during estimation. This can only be used by Windows users for now.
#' @param no_cores Number of cores used if parallel computing chosen. The default is 2. If \code{NULL}, all available cores are used.
#' @param cens.state Code used in the dataset to indicate the censored states.
#' @param living.exact Name of the variable in the dataset indicating whether an observation is exactly observed or not.
#' @param verbose If \code{TRUE}, prints the convergence criterion obtained at each iteration of the full algorithm. The default is \code{FALSE}.
#' @param justComp Can be \code{c('lik', 'grad', 'hess')} (or a subvector of this) to compute the log-likelihood, gradient and Hessian at the starting parameter value, without carrying out fitting. Defaults to \code{NULL}.
#' @param approxHess If \code{TRUE} an approximation of the Hessian based on the outer product of the gradient vector is computed. The default is \code{FALSE}.
#'
#' @importFrom grDevices cm.colors gray heat.colors terrain.colors topo.colors
#' @importFrom graphics lines par text title
#' @importFrom stats pnorm predict printCoefmat qexp qlnorm runif terms uniroot
#' @importFrom utils getFromNamespace
#' @importFrom mgcv Sl.initial.repara Sl.repara exclude.too.far mroot predict.gam totalPenaltySpace
#'
#' @usage fmsm(formula, data, id, state, death, pmethod = 'eigendecomp',
#'        aggregate = TRUE, params.0 = NULL, sp.0 = NULL,
#'        constraint = NULL, sp.method = 'perf', iterlimsp = 50,
#'        Q.diagnostics = TRUE, fit = TRUE, iterlim = 100,
#'        tolsp = 1e-7, tolsp.EFS = 0.1, parallel = FALSE, no_cores = 2,
#'        cens.state = NULL, living.exact = NULL, verbose = FALSE,
#'        justComp = NULL, approxHess = FALSE)
#'
#' @return The function returns an object of class \code{fmsm} as described in \code{fmsmObject}.
#'
#' @examples
#' \dontrun{
#'
#' ##################################################
#' # MULTISTATE SURVIVAL MODELLING with CAV DATA ####
#' ##################################################
#'
#' library(flexmsm)
#'
#' # Import data
#' Data <- IDM_cav
#'
#' # MODEL SPECIFICATION ####
#' formula <- list(years ~ s(years, bs = 'cr', k = 10) + dage + pdiag, # 1-2
#'                 years ~ s(years, bs = 'cr', k = 10) + dage + pdiag, # 1-3
#'                 0,                                                  # 2-1
#'                 years ~ s(years, bs = 'cr', k = 10) + dage + pdiag, # 2-3
#'                 0,                                                  # 3-1
#'                 0                                                   # 3-2
#' )
#'
#' # Counts of pairs of consecutive states observed (C = counts, T = times)
#' counts.CT <- state.pairs.CT(formula = formula, data = Data, time = 'years',
#'                             state = 'state', id = 'PTNUM')
#' counts.CT$counts
#'
#'
#' # MODEL FITTING ###
#'
#' # NOTE ***
#' # Takes about 18 minutes on a machine with Windows 10,
#' # Intel 2.20 GHz core, 16 GB of RAM and 8 cores, using all cores.
#' # The default is to use 2 cores, this takes about 26 minutes.
#' # To use all available cores on your device input no_cores = NULL.
#' # ****
#'
#' fmsm.out <- fmsm(formula = formula, data = Data,
#'                  id = 'PTNUM', state = 'state', death = TRUE,
#'                  fit = TRUE, parallel = TRUE, no_cores = 2,
#'                  pmethod = 'analytic')
#'
#' print(fmsm.out)
#'
#' AIC(fmsm.out)
#' BIC(fmsm.out)
#'
#' # FITTING SUMMARY ####
#' summary(fmsm.out)
#' conv.check(fmsm.out)
#'
#' ####################
#' # VISUALISATION ####
#' ####################
#'
#' # PLOT THE SMOOTHS OF TIME FOR EACH TRANSITION ####
#' # par(mfrow = c(1,3))
#' plot(fmsm.out)
#'
#'
#' # Consider a patient with:
#' dage.pred <- 16      # - 16 year old donor
#' pdiag.pred <- 0      # - IDC as principal diagnosis
#' start.pred <- 0      # - start observation at time t = 0
#' stop.pred <- 15      # - t = 15 years for time horizon
#' n.pred <- 21         # - 21 time points
#' no.state.pred <- -13 # - (because we don't need this, so anything is fine)
#'
#' newdata <- data.frame(PTNUM = rep(1, n.pred),
#'                       years = seq(start.pred, stop.pred, length.out = n.pred),
#'                       state = rep(no.state.pred, n.pred),
#'                       dage = rep(dage.pred, n.pred), pdiag = rep(pdiag.pred, n.pred))
#'
#'
#'
#' # ESTMATED TRANSITION INTENSITIES ####
#'
#' # Plot of estimated transition intensities
#' # par(mfrow = c(1,3))
#' Q.hat <- Q.pred(fmsm.out, newdata = newdata, get.CI = TRUE, plot.Q = TRUE, rug = TRUE,
#'                 ylim = c(0, 1.5))
#'
#' # Estimated transition intensity matrix at, e.g., t = 0
#' round(Q.hat$Q.hist[,,1], 3)
#'
#'
#'
#' # ESTMATED TRANSITION PROBABILITIES ####
#'
#' # Plot of estimated transition probabilities
#' # par(mfrow = c(2,3))
#' P.hat <- P.pred(fmsm.out, newdata = newdata, get.CI = TRUE, plot.P = TRUE, rug = TRUE)
#'
#' # Estimated 15 year transition probability matrix
#' round(P.hat$P.pred, 3)
#' # e.g., there is a 6.2% chance of observing CAV onset 15 years after transplant
#'
#'
#' }
#'
#' @export
#'
#'
fmsm = function(formula, data, id, state, death, pmethod = 'eigendecomp',
                aggregate = TRUE, params.0 = NULL, sp.0 = NULL,
                constraint = NULL, sp.method = 'perf', iterlimsp = 50,
                Q.diagnostics = TRUE, fit = TRUE, iterlim = 100,
                tolsp = 1e-7, tolsp.EFS = 0.1, parallel = FALSE, no_cores = 2,
                cens.state = NULL, living.exact = NULL, verbose = FALSE,
                justComp = NULL, approxHess = FALSE){


  if(parallel == TRUE & Sys.info()['sysname'] != 'Windows') stop('Parallel computing is currently only supported for Windows users. \nPlease contact the authors for further details or change parallel to FALSE.')

  # Added on 04/05/2022 because need this for when we use efs
  Sl.sf = list()


  mod.list = list()

  cov.names = c() # names of covariates included as parametric terms in each transition specification
  cov.n = c() # number of parametric terms in each transition specification
  pcov.mapping = list() # to know which transition has which covariate


  # ******************* #
  # DESIGN MATRIX SETUP #
  # ******************* #

  nstates = 0.5*(1+sqrt(1+4*length(formula))) # this follows from recalling that we need nstates^2 - nstates = off-diag terms = length(formula)

  if( nstates != round(nstates) ) stop('List of model specifications incorrect. You need to provide a formula (or a 0) \n for each non-diagonal transition intensity')

  # Define a matrix indicating the allowed transitions ******
  whereQ = matrix(0, nstates, nstates)
  k = 1
  l = 1
  for(i in 1:nstates){
    for(j in 1:nstates){
      if(i != j) {
        if( formula[[k]] != 0) {whereQ[i,j] = l; l = l+1} # instead of just 1, put progressive number, useful for later (when setting up full Q matrix)
        k = k+1
      }
    }
  }
  rm(k)
  rm(l)

  # warning 1: analytic only available for IDMs for now
  if(!(sum(whereQ != 0) == 3 & whereQ[1,2] != 0 & whereQ[1,3] != 0 & whereQ[2,3] != 0) & pmethod == 'analytic'){
    warning('pmethod cannot be \'analytic\', only IDM process is supported for now.\n pmethod has been changed to \'eigendecomp\'.')
    pmethod = 'eigendecomp'
  }
  # ***********************************************************

  # cl <- match.call() # not sure if this is needed??
  full.mf <- match.call(expand.dots = FALSE)

  # # old code
  # if(is.null(living.exact)){
  #   ind = match(c('data', 'id', 'state'), names(full.mf), nomatch = 0)
  # }  else {
  #   ind = match(c('data', 'id', 'state', 'living.exact'), names(full.mf), nomatch = 0)
  # }
  #
  # mf = full.mf[c(1, ind)] # mf drops all other arguments (only for the purpose of obtaining dataset)

  # new code ****
  ind.df = match(c('data'), names(full.mf), nomatch = 0)
  # ind.id = match(c('id'), names(full.mf), nomatch = 0)
  # ind.st = match(c('state'), names(full.mf), nomatch = 0)
  #
  # if(!is.null(living.exact)) ind.LE = match(c('living.exact'), names(full.mf), nomatch = 0)

  mf = full.mf[c(1,ind.df)] # mf drops all other arguments (only for the purpose of obtaining dataset)
  # ****


  short.formula = formula
  i = 1
  while(i <= length(short.formula)) {if(short.formula[[i]] == 0) {short.formula[[i]] = NULL} else{i = i + 1} }
  rm(i)

  l.short.formula = length(short.formula)



  ig = mgcv::interpret.gam(short.formula)

  # *** Define the internal time-to-event ***
  tte = ig$response
  # *****************************************

  fake.formula = ig$fake.formula

  mf$formula = fake.formula
  mf[[1]] = as.name('model.frame')

  # mf$drop.unused.levels = TRUE # perhaps useful for neater objects ?

  og.data = data
  data = eval(mf, parent.frame())

  # new code ****
  data$'(id)' = og.data[,id] #[,full.mf[ind.id][[1]]]
  data$'(state)' = og.data[,state] #[,full.mf[ind.st][[1]]]
  if(!is.null(living.exact)) data$'(livingexact)' = og.data[, living.exact] #full.mf[ind.LE][[1]]]
  # ****

  # ********************************
  # Now it's time to transform the data in the form we want it to be for LikGradHess
  # Note sure if to add this to the data itself or to full.X. I think the former is more
  # appropriate because in the latter we have the actual design matrices associated to each
  # transition while the info in which transition occurred is a "global" one... think about this.
  # Note2: eventually to aggregation as well
  # Note3: modify original dataset or mantain the old one and create a copy? For now MODIFY EXISTING.
  # Note4: need to omit data here because when I evaluate the Q matrices on the data I only need it
  # in the left extremity of every interval for every individual (so last time should never be used).

  data$"(tostate)" = data$"(fromstate)" = NA # note: in the following the NA is put AFTER to follow the logic that in piece-wise const approx we take the value from the left extremity of the interval

  # !!! NOTE FOR POTENTIAL IMPROVEMENT: PERHPAPS TO THE FOLLOWING USING duplicate() ON THE ID... THINK ABOUT THIS !!!
  for(id.t in unique(data$"(id)")){
    data$"(fromstate)"[data$"(id)" == id.t] = c(data$"(state)"[data$"(id)" == id.t][-sum(data$"(id)" == id.t)], NA) # remove last one
    data$"(tostate)"[data$"(id)" == id.t] = c(data$"(state)"[data$"(id)" == id.t][-1], NA) # remove first one
    data$"(timelag)"[data$"(id)" == id.t] = c(diff(data[data$"(id)" == id.t, tte]), NA)
  }

  if (is.null(data$"(living.exact)")) data$"(living.exact)" = FALSE

  n = nrow(data)
  N = length(unique(data$"(id)"))


  # # **************
  # data.long = data
  # # **************
  # # data = na.omit(data) # am I removing any necessary information like this?
  # # ********************************




  # CREATE AGGREGATED DATASET ************************
  # Two or more rows can be aggregated when they have the same values fpr
  # - covariates and timelag (imply same Q matrix and P matrix)
  # - fromstate and tostate (imply, in particular, same likelihood contribution)

  # Get covariate names
  only.cov = ig$pred.formula
  only.cov = attr(terms(only.cov), 'term.labels')

  # In the following there are no duplicate rows... BUT NEED WAY TO COUNT DUPLICATE ROWS
  data.agg = data[!duplicated(data[, c('(fromstate)', '(tostate)', '(timelag)', only.cov) ], MARGIN = 2),]
  pasted.data =  do.call('paste', data[, c('(fromstate)', '(tostate)', '(timelag)', only.cov)])
  counts = as.numeric(table(pasted.data))

  # Create new variable in which the number of repetitions for given row will be saved
  data.agg$'(nrep)' = NA
  data.agg$'(nrep)'[order(unique(pasted.data))] = counts


  # Not sure if this is needed, but perhaps add (nrep) also to non aggregated datasets? (just set it all
  # to 1). Also, more brutal aggregation can take place if we drop fromstate and tostate (indeed the
  # problem lies mainly in having to compute the full transition probability matrix, which changes only as
  # the covariates and timelag change, i.e. we compute it for all fromstate and tostate regardless whether it
  # is repeated or not). THINK ABOUT THIS.

  data$'(nrep)' = 1

  # *****************************************************



  # *****************************************************************
  # *****************************************************************
  # Use gam to build design matrix for each transition intensity
  # (not efficient but only thing I can think of until I understand
  # how to unite formulas in one global formula).
  # *****************************************************************
  # *****************************************************************

  # *** Long version w/ intercept *** for comparison with msm (i.e. do not remove last observation for each individual)

  full.X = c()

  start.pos.par = c()
  start.pos.par.only.smooth = c()
  start.pos.par.only.smooth.FPC = c() # FPC = for penalty construction
  start.pos.par.detailed = c()

  S.list = list()
  sp.counter = 1 # NEW [30/03/2022]

  n.smooth = c()

  for(i in 1:length(short.formula)){
    mod = mgcv::gam(short.formula[[i]], data = data) # NOTE! Here we are using long data, this is just per far tornare le cose con la costruzione degli smooth done in old benchmark code (and also because this is what we do in general)

    # save the mods just for aiding prediction temporarily - [16/05/2022] ****
    assign(paste('mod', i, sep = ''), mod)
    # ************************************************************************

    # [04/05/2022] This is needed for efs smoothing estimation method ******************
    mod.for.setup = mgcv::gam(short.formula[[i]], data = data, fit = FALSE)
    Sl.sf[[i]] = mgcv::Sl.setup(mod.for.setup)
    # **********************************************************************************

    matr = predict(mod, type = 'lpmatrix')
    assign(paste('X', i, sep = ''), matr) # WITH INTERCEPT - [comment from 30/03/2022] this was originally here just in case we wanted to keep the design matrices relating to each transition intensity separate, do we want to keep this?


    # For penalty setup
    # *****************
    if(length(mod$smooth) > 0){ # 30/03/2022 - change added now that we are transitioning to more complex models, this will be the new norm but keeping it separate now just for safety
      for(ii in 1:length(mod$smooth)){
        assign(paste('pen', ii, sep = ''),  mod$smooth[[ii]]$S[[1]])  # retain just the first pen per smooth - anyway won't be using full.S anymore for construction of penalty, just gives dimensions now
        for(ii.inner in 1:length(mod$smooth[[ii]]$S)){
           # do we ever have S[[j]] with j != 1 ?? YES! When we  have tensor interactions, then we have one S per directions and one sp per S as well!!
          S.list[[sp.counter]] = mod$smooth[[ii]]$S[[ii.inner]]
          sp.counter = sp.counter + 1
        }
      }
    }
    # ******************



    n.smooth = c(n.smooth, length(mod$smooth))
    pterms.temp = attr(mod$pterms, 'term.labels')
    n.pterms = length(pterms.temp)

    # Saving linear covariate names so constraints can be put here ****
    cov.names = c(cov.names, pterms.temp)
    cov.n = c(cov.n, n.pterms)

    if(n.pterms > 0){
      # To add the positions of covariates found
      for(i.cov in 1:n.pterms){
        if( pterms.temp[i.cov] %in% names(pcov.mapping) ){
          pcov.mapping[[pterms.temp[i.cov]]] = c(pcov.mapping[[pterms.temp[i.cov]]], 1)
        } else {
          pcov.mapping[[pterms.temp[i.cov]]] = c(rep(0,i-1), 1)
        }
      }
    }

    # To register when a covariate is not present for the given transition
    if(length(names(pcov.mapping)) > 0 & any(!(names(pcov.mapping) %in% pterms.temp))){
      not.pres.idx = which(!(names(pcov.mapping) %in% pterms.temp))
      for(i.map in 1:length(not.pres.idx)){
        pcov.mapping[[not.pres.idx[i.map]]] = c(pcov.mapping[[not.pres.idx[i.map]]], 0)
      }
    }


    # *****************************************************************

    if(i == 1) { # record position of parameters with respect to global design matrix
      start.pos.par = c(start.pos.par, 1, ncol(matr)+1)
      start.pos.par.detailed = c(start.pos.par.detailed, 1, n.pterms+2)

      # For penalty setup ***
      if (length(mod$smooth) > 0){
        pen = get(paste('pen', 1, sep = ''))
        full.S = cbind( matrix(0, ncol = n.pterms, nrow = nrow(pen)+n.pterms),
                             rbind(matrix(0, ncol = ncol(pen), nrow = n.pterms), pen) ) # add zeroes matrix for parametric terms

        start.pos.par.detailed = c(start.pos.par.detailed, ncol(full.S)+2) # plus 2 because we add intercept before everything but this happens later on in code
        start.pos.par.only.smooth = c(start.pos.par.only.smooth, n.pterms+2)

        if(length(mod$smooth[[1]]$S) == 1) {
          start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, n.pterms+2)
        } else if(length(mod$smooth[[1]]$S) == 2) {
          start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, n.pterms+2, n.pterms+2) # repeat twice --- can we have more than two?? can there be three way tensor interactions?? For now no...
        } else {
          stop('There seem to be three penalty matrices associated with a single spline. Please contact authors for further details.')
        }


        if(length(mod$smooth) > 1){ # 16/09/2022 - because this will be the general code for when we also just have the smooth of time and nothing else
          for(ii in 2:length(mod$smooth)){
            start.pos.par.only.smooth = c(start.pos.par.only.smooth, ncol(full.S)+2) # CHECK IF HOLDS TRUE ALSO FOR MORE SPLINES

            if(length(mod$smooth[[ii]]$S) == 1) {
              start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, ncol(full.S)+2)
            } else if(length(mod$smooth[[ii]]$S) == 2) {
              start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, ncol(full.S)+2, ncol(full.S)+2) # repeat twice --- can we have more than two?? can there be three way tensor interactions?? For now no...
            } else {
              stop('There seem to be three penalty matrices associated with a single spline. Please contact authors for further details.')
            }

            pen = get(paste('pen', ii, sep = ''))
            full.S = cbind( rbind(full.S, matrix(0, ncol = ncol(full.S), nrow = nrow(pen))),
                                 rbind(matrix(0, ncol = ncol(pen), nrow = nrow(full.S)), pen)) # add zeroes matrix for parametric terms
            start.pos.par.detailed = c(start.pos.par.detailed, ncol(full.S)+2) # (same here) # plus 2 because we add intercept before everything but this happens later on in code
          }
        }


      } else {
        full.S = matrix(0, ncol = n.pterms, nrow = n.pterms)
        # start.pos.par.detailed = c(start.pos.par.detailed, ncol(full.S)+1)

      }

      full.S = cbind( rep(0, nrow(full.S)+1), rbind(rep(0, ncol(full.S)), full.S) )  # add zeroes for intercept

      # ************************

    } else { # otherwise we record also the position of the very last parameter
      start.pos.par = c(start.pos.par, start.pos.par[length(start.pos.par)] + ncol(matr) )

      # For penalty setup ***
      full.S = cbind( rbind(full.S, rep(0, ncol(full.S))), rep(0, nrow(full.S)+1) )  # add zeroes for intercept
      if (length(mod$smooth) > 0){
        full.S = cbind( rbind(full.S, matrix(0, ncol = ncol(full.S), nrow = n.pterms)),
                             matrix(0, ncol = n.pterms, nrow = nrow(full.S)+n.pterms) ) # add zeroes matrix for parametric terms

        start.pos.par.detailed = c(start.pos.par.detailed, ncol(full.S)+1)
        start.pos.par.only.smooth = c(start.pos.par.only.smooth, ncol(full.S)+1)

        if(length(mod$smooth[[1]]$S) == 1) {
          start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, ncol(full.S)+1)
        } else if(length(mod$smooth[[1]]$S) == 2) {
          start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, ncol(full.S)+1, ncol(full.S)+1) # repeat twice --- can we have more than two?? can there be three way tensor interactions?? For now no...
        } else {
          stop('There seem to be three penalty matrices associated with a single spline. Please contact authors for further details.')
        }


        for(ii in 1:length(mod$smooth)){

          pen = get(paste('pen', ii, sep = ''))
          full.S = cbind( rbind(full.S, matrix(0, ncol = ncol(full.S), nrow = nrow(pen))),
                               rbind(matrix(0, ncol = ncol(pen), nrow = nrow(full.S)), pen))

          start.pos.par.detailed = c(start.pos.par.detailed, ncol(full.S)+1)

          if(ii < length(mod$smooth)){

            start.pos.par.only.smooth = c(start.pos.par.only.smooth, ncol(full.S)+1)

            if(length(mod$smooth[[ii+1]]$S) == 1) { # has to be +1 because at each ii we are adding the start pos of the *next* smooth so if its double add twice
              start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, ncol(full.S)+1)
            } else if(length(mod$smooth[[ii+1]]$S) == 2) {
              start.pos.par.only.smooth.FPC = c(start.pos.par.only.smooth.FPC, ncol(full.S)+1, ncol(full.S)+1) # repeat twice --- can we have more than two?? can there be three way tensor interactions?? For now no...
            } else {
              stop('There seem to be three penalty matrices associated with a single spline. Please contact authors for further details.')
            }
          }

          if(ii == length(mod$smooth) & i == length(short.formula)) start.pos.par.only.smooth = c(start.pos.par.only.smooth, ncol(full.S)+1) # this is to add the termination point+1 of the smooths
        }

      } else {
        full.S = cbind( rbind(full.S, matrix(0, ncol = ncol(full.S), nrow = n.pterms)),
                             matrix(0, ncol = n.pterms, nrow = nrow(full.S)+n.pterms) ) # add zeroes matrix for parametric terms

        start.pos.par.detailed = c(start.pos.par.detailed, ncol(full.S)+1)
      }


      # *********************
    }

    full.X = cbind( full.X, matr ) # WITH INTERCEPT

    mod.list[[i]] = mod

  }

  # Drop repeated columns (might be useful)
  unique.full.X = unique(full.X, MARGIN = 2)



  # *** Long AGGREGATED version w/ intercept ***

  full.X.agg = full.X[!duplicated(data[, c('(fromstate)', '(tostate)', '(timelag)', only.cov)]), ]

  start.pos.par.agg = start.pos.par
  start.pos.par.only.smooth.agg = start.pos.par.only.smooth
  S.list.agg = S.list

  # Drop repeated columns (might be useful)
  unique.full.X.agg = unique(full.X.agg, MARGIN = 2)

  # ***********************************************************************************************************


  # warning 2 - the formula inputted is not coherent with the data
  # tmp.time = formula[[1]][[2]]  # this is just to get the time variable name (not as a string)
  # tmp.state = full.mf[match('state', names(full.mf), nomatch = 0)][[1]]
  # tmp.id = full.mf[match('id', names(full.mf), nomatch = 0)][[1]]

  idx.lab = which(t(whereQ) != 0, arr.ind = T)

  chck.counts = statetable.fmsm(state = state, subject = id, data = og.data)
  if(nrow(chck.counts) < ncol(chck.counts)) chck.counts = rbind(chck.counts, rep(0, ncol(chck.counts)))
  chck.countsVec = t(chck.counts)[t(whereQ) != 0]

  if(any(chck.countsVec == 0)){
    wrong.trans = which(chck.countsVec == 0)
    for(wr.tr in wrong.trans){
      if(any(chck.counts[1:idx.lab[wr.tr, 2], idx.lab[wr.tr, 1]] != 0) & is.null(params.0)) stop('An equation for transition ', idx.lab[wr.tr, 2], '->', idx.lab[wr.tr, 1],
                                                                              ' was specified, this direct transition \nis never observed in the data but a path is possible. Please input starting parameters manually through params.0 argument.')
      if(all(chck.counts[1:idx.lab[wr.tr, 2], idx.lab[wr.tr, 1]] == 0)) stop('An equation for transition ', idx.lab[wr.tr, 2], '->', idx.lab[wr.tr, 1],
                                                                             ' was specified, but this path \nis never observed in the data. Please replace with a 0 in the formula.')
    }
  }

  # Retrieve starting values for parameters
  if(is.null(params.0)){
    # initQ = eval(substitute(msm::crudeinits.msm(formula = tmp.state ~ tmp.time, subject = tmp.id, qmatrix = whereQ, data = og.data),
    #                         list(tmp.state = tmp.state, tmp.time = tmp.time, tmp.id = tmp.id)))
    initQ = crudeinits.fmsm(state = state, time = tte, subject = id, qmatrix = whereQ, data = og.data)
    initQ = log(initQ[whereQ != 0])

    params.0 = c()
    for(i.p0 in 1:length(short.formula)) params.0 = c(params.0, initQ[i.p0], rep(0, ncol(eval(parse(text = paste('X', i.p0, sep = ''))))-1))
  }

  # Retrieve starting values for smoothing parameters
  if(is.null(sp.0)) sp.0 = exp(rep(2, length(S.list)))



  # MAKING GENERAL CODE FOR CONSTRAINING PART ***************
  if(is.null(constraint)) constraint = list()

  covlabels = unique(cov.names)
  idx.pcov = c(0, cumsum(cov.n))

  # # Complicated P-splines I model from Mariano thesis
  # constraint = list(time2 = c(1, 2, 3, 2, 3, 2, 3, 2),
  #                   sex = c(1, 2, 3, 2, 4, 2, 2))


  # constraint = list(dage = c(1, 1, 1),
  #                   pdiag = c(1, 1, 1))


  if(length(covlabels) > 0){
    for(i in 1:length(covlabels)){
      if( !(covlabels[i] %in% names(constraint))  ){
        constraint[[covlabels[i]]] = 1:max(whereQ)
      }

      if( length(constraint[[covlabels[i]]]) < max(whereQ)  ){
        constr.temp = rep(1, max(whereQ))
        constr.temp[which(pcov.mapping[[covlabels[i]]] == 1)] = constraint[[covlabels[i]]]+1
        constraint[[covlabels[i]]] = constr.temp
        rm(constr.temp)
      }
    }

    all.terms = attr(full.X, 'dimnames')[[2]]
    for(q.term in 1:max(whereQ)){

      if(q.term == 1){
        pos.optparams2 = pos.optparams = start.pos.par[q.term]:(start.pos.par[q.term+1]-1)
        count2 = count = pos.optparams[length(pos.optparams)]
      } else {

        q.term.names = all.terms[start.pos.par[q.term]:(start.pos.par[q.term+1]-1)]

        for(q.term.single in q.term.names){

          if( q.term.single %in% names(constraint) ){ # so essentially if it is a parametric term
            if( any(constraint[[q.term.single]][q.term] %in% constraint[[q.term.single]][1:(q.term-1)]) ){ #then this is constrained to an existing thing

              # trace back the thing through the names and positions already recorded
              which.base.q = which(constraint[[q.term.single]][1:(q.term-1)] == constraint[[q.term.single]][q.term])[1] # this tells me which transition it is in (the base param)
              constr.pos = pos.optparams[start.pos.par[which.base.q] +  match(q.term.single, all.terms[start.pos.par[which.base.q]:(start.pos.par[which.base.q+1]-1)]) - 1]

              pos.optparams = c(pos.optparams, constr.pos)
              pos.optparams2 = c(pos.optparams2, constr.pos)

              # count = count + 0 # do not update this one
              count2 = count2 + 1 # update

            } else {
              pos.optparams = c(pos.optparams, count + 1)
              count = count + 1

              pos.optparams2 = c(pos.optparams2, count2 + 1)
              count2 = count2 + 1
            }
          } else { # so it is a smooth term (since all parametric terms are recorded at this point in the constraints list) or intercept
            pos.optparams = c(pos.optparams, count + 1)
            count = count + 1

            pos.optparams2 = c(pos.optparams2, count2 + 1)
            count2 = count2 + 1
          }

        }
      }
    }

  } else {
    pos.optparams = pos.optparams2 = 1:length(params.0)
  }




  params.0 = params.0[1:length(params.0) == pos.optparams2] # this will return the same params.0 if there are no constraints


  # *********************************************************

  if(length(sp.0) != (sp.counter-1) ) stop(paste('The length smoothing parameters vector provided (', length(sp.0),') does not correspond with the number of smooths found (', sp.counter-1, ').'), sep = '')
  if(length(params.0) != max(pos.optparams) ) stop(paste('The length parameters vector provided (', length(params.0),') does not correspond with the model specification provided. Should contain', max(pos.optparams), 'elements.'), sep = '')



  # Gathering the necessary variable to be fed into model fitting
  suStf = list(data = data, nstates = nstates,
               start.pos.par = start.pos.par, pos.optparams = pos.optparams,
               pos.optparams2 = pos.optparams2, start.pos.par.only.smooth = start.pos.par.only.smooth,
               start.pos.par.only.smooth.FPC = start.pos.par.only.smooth.FPC,
               start.pos.par.only.smooth.FPC.constr = pos.optparams[start.pos.par.only.smooth.FPC],
               start.pos.par.only.smooth.constr = pos.optparams[start.pos.par.only.smooth], # need this for sp estimation when constraints in place, otherwise coincides with start.pos.part.only.smooth. NOTE: the last element is always NA here just because by convention the OG param has as last element (its length+1), which will always not exist when selecting with pos.optparams
               start.pos.par.constr = c(pos.optparams[start.pos.par[-length(start.pos.par)]], max(pos.optparams)+1), # need this as more general version of start.pos.par, coincides with it when no constraining present - perhaps just keep this in the future
               start.pos.par.detailed = start.pos.par.detailed,
               l.short.formula = l.short.formula,
               whereQ = whereQ, full.X = full.X, Sl.sf = Sl.sf, S.list = S.list, full.S = full.S,
               sp = sp.0, params = params.0,
               cens.state = cens.state,
               tte = tte,
               death = death, pmethod = pmethod,
               constraint = constraint,
               mf = mf, mod.list = mod.list)




  msm.fit.object = msm.post.object = logLik = t.edf = NULL
  singleComp = NULL


  if(!is.null(justComp)){ # computation just of a single lik,grad,hess (at starting parameter)

    if(fit) stop('justComp can only be used when fit = FALSE, i.e. the model fitting need not be carried out and \nonly one instance of the likelihood (gradient and/or hessian) needs to be computed at the \nparameters provided in params.0.')

    singleComp = msmJustOne(params.0 = params.0, sp = sp.0, pmethod = pmethod, suStf = suStf, death = death,
                            Q.diagnostics = Q.diagnostics, parallel = parallel, no_cores = no_cores, justComp = justComp)

  }



  if(fit){

    # **************** #
    # MODEL FITTING ####
    # **************** #
    msm.fit.object = msm.fit(params.0 = params.0, sp = sp.0, pmethod = pmethod,
                             suStf = suStf, death = death,
                             Q.diagnostics = Q.diagnostics, iterlimsp = iterlimsp,
                             iterlim = iterlim,
                             sp.method = sp.method, approxHess = approxHess,
                             verbose = verbose, tolsp = tolsp, tolsp.EFS = tolsp.EFS, parallel = parallel, no_cores = no_cores)
    # ********************** #
    # POST FITTING THINGS ####
    # ********************** #
    msm.post.object = msm.fit.post(msm.fit.object = msm.fit.object, mod.list = mod.list, suStf = suStf)

    logLik = msm.post.object$logLik
    t.edf = msm.post.object$t.edf

  }



  L = list(suStf = suStf,
       # unique.full = unique.full.X, # do we even need this?
       msm.fit.object = msm.fit.object,
       msm.post.object = msm.post.object,
       formula = formula,
       short.formula = short.formula,
       n = n, N = N,
       logLik = logLik, t.edf = t.edf,
       singleComp = singleComp
  )
  class(L) <- c("fmsm")
  L

}

