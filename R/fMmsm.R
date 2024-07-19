#' Flexible transition intensity based models for two dependant multistate processes
#'
#' @description
#' XXXX.
#'
#' @param formula1 Model specification for the transition intensities of the first process.
#' @param data1 Dataset of the first process.
#' @param params1 XXX.
#' @param spP1 Smoothing parameter for the first process.
#' @param constraint1 XXX.
#' @param id1 Name of the variable in the dataset representing the unique code associated with each patient in the first process.
#' @param state1 Name of the variable in the first process dataset representing the state occupied by the patient at the given time.
#' @param formula2 Model specification for the transition intensities of the second process.
#' @param data2 Dataset of the second process.
#' @param params2 XXX.
#' @param spP2 Smoothing parameter for the second process.
#' @param constraint2 XXX.
#' @param id2 Name of the variable in the dataset representing the unique code associated with each patient in the second process.
#' @param state2 Name of the variable in the second process dataset representing the state occupied by the patient at the given time.
#' @param phi XXX.
#' @param parallel If \code{TRUE} parallel computing is used during estimation. This can only be used by Windows users for now.
#' @param no_cores Number of cores used if parallel computing chosen. The default is 2. If \code{NULL}, all available cores are used.
#' @param iterlim Maximum allowed iterations for trust region algorithm.
#' @param verbose XXX.
#' @param pmethod Which method should be used for the computation of the transition probability matrix. Available options are
#' \itemize{
#' \item \code{'eigendecomp'} (default): this method is based on the eigendecomposition of the transition intensity matrix (from Kalbfleisch & Lawless 1985);
#' \item \code{'analytic'}: uses analytic expressions of the transition probabilities, obtained by solving the Kolmogorov forward differential equation, only implemented for IDMs for now;
#' \item \code{'scaling&squaring'}: this is the scaling and squaring method implemented as proposed in Fung (2004).This is inefficient, so its use is not recommended. Can be used to investigate convergence errors.
#' }
#' @param aggregate Whether or not data should be aggregated (this slightly improves efficiency as redundancies in the data are eliminated). The default is \code{TRUE}.
#' @param sp.method Method to be used for smoothing parameter estimation. The default is \code{magic}, the automatic multiple smoothing parameter selection algorithm. Alternatively, \code{efs} can be used for the Fellner-Schall method. To suppress the smoothing parameter estimation set this to \code{NULL}.
#' @param iterlimsp Maximum allowed iterations for smoothing parameter estimation.
#' @param Q.diagnostics If \code{TRUE}, diagnostics information on the Q matrix are saved. The default \code{TRUE}.
#' @param tolsp Convergence criterion used in \code{magic} based smoothing parameter estimation.
#' @param tolsp.EFS Convergence criterion used in \code{efs} based smoothing parameter estimation.

#'
#'
#' @usage fMmsm(formula1, data1, id1, state1,
#'        params1 = NULL, spP1 = NULL, constraint1 = NULL,
#'        formula2, data2, id2, state2,
#'        params2 = NULL, spP2 = NULL, constraint2 = NULL,
#'        phi = NULL,
#'        pmethod = 'eigendecomp',
#'        aggregate = TRUE, sp.method = 'perf', iterlimsp = 50,
#'        Q.diagnostics = TRUE, iterlim = 100, verbose,
#'        tolsp = 1e-7, tolsp.EFS = 0.1, parallel = FALSE, no_cores = 2)
#'
#' @return The function returns an object of class \code{fmsm} as described in \code{fmsmObject}.
#'
#'
#' @export
#'
#'
#'
fMmsm = function(formula1, data1, id1, state1,
                 params1 = NULL, spP1 = NULL, constraint1 = NULL,
                 formula2, data2, id2, state2,
                 params2 = NULL, spP2 = NULL, constraint2 = NULL,
                 phi = NULL,
                 pmethod = 'eigendecomp',
                 aggregate = TRUE, sp.method = 'perf', iterlimsp = 50,
                 Q.diagnostics = TRUE, iterlim = 100, verbose,
                 tolsp = 1e-7, tolsp.EFS = 0.1, parallel = FALSE, no_cores = 2){




  # (TO DO: include a check that the two processes have the same number of individuals)


  # Starting values for the parameters of PROCESS 1
  if(is.null(params1) | is.null(spP1)){
    prelim1 <- fmsm(formula = formula1, data = data1,
                   id = id1, state = state1, death = FALSE,
                   fit = FALSE, justComp = NULL)
    if(is.null(params1)) params1 <- prelim1$suStf$params
    if(is.null(spP1)) spP1 <- prelim1$suStf$sp
  }

  # Starting values for the parameters of PROCESS 2
  if(is.null(params2) | is.null(spP1)){
    prelim2 <- fmsm(formula = formula2, data = data2,
                   id = id2, state = state2, death = FALSE,
                   fit = FALSE, justComp = NULL)
    if(is.null(params2)) params2 <- prelim2$suStf$params
    if(is.null(spP2)) spP2 <- prelim2$suStf$sp
  }

  # Starting values for the parameters of the COPULA
  if(is.null(phi)) phi <- 2 # small-ish number if the user does not provide the copula parameter

  paramsMulti <- c(params1, params2, phi)




  # **************** #
  # MODEL FITTING ####
  # **************** #
  Mmsm.fit.object = Mmsm.fit(paramsMulti = paramsMulti,
                             formula1 = formula1, data1 = data1,
                             id1 = id1, state1 = state1, spP1 = spP1, l.params1 = length(params1),
                             formula2 = formula2, data2 = data2,
                             id2 = id2, state2 = state2,  spP2 = spP2, l.params2 = length(params2),
                             pmethod = pmethod,
                             Q.diagnostics = Q.diagnostics, iterlimsp = iterlimsp,
                             iterlim = iterlim,
                             sp.method = sp.method,
                             verbose = verbose, tolsp = tolsp, tolsp.EFS = tolsp.EFS,
                             parallel = parallel, no_cores = no_cores)

  # ********************** #
  # POST FITTING THINGS ####
  # ********************** #
  Mmsm.post.object = Mmsm.fit.post(Mmsm.fit.object = Mmsm.fit.object, l.params1 = length(params1), l.params2 = length(params2))

  logLik = Mmsm.post.object$logLik
  t.edf = Mmsm.post.object$t.edf



  L <- list(Mmsm.fit.object = Mmsm.fit.object,
            Mmsm.post.object = Mmsm.post.object)


  class(L) <- c("fMmsm")
  L



}
