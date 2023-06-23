
#' @title Predict and plot the transition intensities
#'
#' @description
#' Function to predict and plot the estimated transition intensities (and confidence intervals).
#'
#' @param object Fitted model object.
#' @param newdata Dataframe containing the profile for which one wished to obtain the predicted transition intensities.
#' @param get.CI Whether to compute the confidence intervals.
#' @param n.sim.CI Number of simulations to be used for confidence intervals computation.
#' @param prob.lev Probability level of confidence intervals.
#' @param plot.Q Whether to output plots of transition intensities.
#' @param which.plots Number between 1 and the maximum number of non-null transition intensities. This can be used if only some plots are to be plotted.
#' @param cond.list.2d Value of covariate(s) to be kept fixed in the plotting of 3D-transition intensities.
#' @param plot.Q.2d Whether to plot 3D transition intensities (only valid if 2D-smooths are present).
#' @param rug Whether to include a rugplot of the observed transition times.
#' @param ... Other graphical parameters.
#'
#' @usage Q.pred(object, newdata, get.CI = TRUE,
#'        n.sim.CI = 1000, prob.lev = 0.05,
#'        plot.Q = FALSE, which.plots = NULL,
#'        cond.list.2d = NULL, plot.Q.2d = FALSE,
#'        rug = TRUE, ...)
#'
#' @return Estimated transition intensities (and confidence intervals).
#' \item{Q.hist}{List of predicted transition intensity matrices computed at each time point specified in \code{newdata}. This is a \code{nstates x nstates x n.pred} array, where \code{n.pred} is the number of rows in \code{newdata}.}
#' \item{Q.CI.lower}{Matrix containing the lower bounds of the confidence intervals for the predicted transition intensity matrix.}
#' \item{Q.CI.upper}{Matrix containing the upper bounds of the confidence intervals for the predicted transition intensity matrix.}
#' \item{full.X}{Full design matrix corresponding to the \code{newdata} provided.}
#' \item{Q.sim.hist}{List of transition intensity matrices simulated to obtain the confidence intervals at each time point from \code{newdata}. May be useful to quickly obtain intervals for a different confidence level.}
#'
#'
#' @export
#'
#' @seealso \code{\link{fmsm}}
#'

Q.pred = function(object, newdata,
                      get.CI = TRUE, n.sim.CI = 1000,
                      prob.lev = 0.05, plot.Q = FALSE,
                      which.plots = NULL, cond.list.2d = NULL, plot.Q.2d = FALSE, rug = TRUE, ...){



  Q.CI.hist = Q.CI.lower = Q.CI.upper = NULL

  params.hat = object$msm.fit.object$fit$argument
  sp.hat = object$msm.fit.object$sp

  n.pred = nrow(newdata)
  if(n.pred < 2) stop('newdata needs to have at least two rows.')

  suStf = object$suStf
  pmethod = suStf$pmethod
  death = suStf$death

  smooth.2d.idx = list()

  full.X = c()
  for(i in 1:length(object$msm.post.object$msm.post.object)){
    mod = object$msm.post.object$msm.post.object[[i]]
    collecting = c()

    for(ii.2d in 1:length(mod$smooth)){
      if(length(mod$smooth[[ii.2d]]$term) == 2){
        collecting = c(collecting, ii.2d)
      }  else if(length(mod$smooth[[ii.2d]]$term) > 2) {
        warning('Unsupported smooth type for plotting. Please contact the authors.')
      }
    }

    smooth.2d.idx[[i]] = collecting

    matr = predict(object = mod, newdata = newdata, type = 'lpmatrix')
    full.X = cbind( full.X, matr)
  }

  # Setup Q matrix
  Qmatr = Q.matr.setup.general(params = params.hat[suStf$pos.optparams], nstates = suStf$nstates,
                               full.X = full.X, start.pos.par = suStf$start.pos.par,
                               l.short.formula = suStf$l.short.formula, whereQ = suStf$whereQ,
                               firstD = FALSE, secondD = FALSE, bound.eta = FALSE,
                               pos.optparams = suStf$pos.optparams, pos.optparams2 = suStf$pos.optparams2)$Qmatr

  Qmatr.debug = Qmatr


  idx.lab = which(t(suStf$whereQ) !=0, arr.ind = T)


  if(get.CI){

    # Confidence intervals ####
    beta.mu = params.hat
    beta.sigma = matrix(object$msm.post.object$Vb, nrow = length(beta.mu), ncol = length(beta.mu))
    # set.seed(24)
    bs = GJRM::rMVN(n = n.sim.CI, mean = beta.mu, sigma = beta.sigma)

    Q.CI.hist = array(dim = c(suStf$nstates, suStf$nstates, n.pred, n.sim.CI))

    for(i.CI in 1:n.sim.CI){

      Q.CI.hist[, , , i.CI] = Q.matr.setup.general(params = bs[i.CI,][suStf$pos.optparams], nstates = suStf$nstates,
                                   full.X = full.X, start.pos.par = suStf$start.pos.par,
                                   l.short.formula = suStf$l.short.formula, whereQ = suStf$whereQ,
                                   firstD = FALSE, secondD = FALSE, bound.eta = FALSE,
                                   pos.optparams = suStf$pos.optparams, pos.optparams2 = suStf$pos.optparams2)$Qmatr
    }


    Q.CI.lower = Q.CI.upper = array(0, dim = c(suStf$nstates, suStf$nstates, n.pred))

    for(ii in 1:sum(suStf$whereQ != 0)){
      Q.CI.lower[idx.lab[ii, 2], idx.lab[ii, 1],] = matrixStats::rowQuantiles(Q.CI.hist[idx.lab[ii, 2], idx.lab[ii, 1],,], probs = c(prob.lev/2), na.rm = TRUE)
      Q.CI.upper[idx.lab[ii, 2], idx.lab[ii, 1],] = matrixStats::rowQuantiles(Q.CI.hist[idx.lab[ii, 2], idx.lab[ii, 1],,], probs = c(1 - prob.lev/2), na.rm = TRUE)
    }

    # if(i.CI/n.sim.CI*100 %% 20 == 0) print(paste(i.CI/n.sim.CI*100, '% of simulations for the CIs completed.'))

  }


  if(plot.Q){

    if(is.na(match('xlab', names(match.call())))) xlab = suStf$tte
    if(is.na(match('lwd', names(match.call())))) lwd = 1
    if(is.na(match('cex.lab', names(match.call())))) cex.lab = 1.4
    if(is.na(match('cex.axis', names(match.call())))) cex.axis = 1.25

    times = newdata[, suStf$tte]
    which.plots = ifelse(is.null(which.plots), sum(suStf$whereQ != 0), which.plots)

    for(ii in 1:which.plots){

      if(is.na(match('ylim', names(match.call())))) ylim = c(0, max(Q.CI.upper[idx.lab[ii, 2], idx.lab[ii, 1],]))


      plot(times, Qmatr[idx.lab[ii, 2], idx.lab[ii, 1],], type = 'l', lwd = lwd,
           xlab = '',
           ylab = paste(idx.lab[ii, 2], 'to', idx.lab[ii, 1], 'trans. intensity'),
           cex.lab = cex.lab, cex.axis = cex.axis, ...)
      title(xlab = xlab, line = 4, cex.lab = cex.lab)


      if(get.CI) {
        lines(times, Q.CI.lower[idx.lab[ii, 2], idx.lab[ii, 1],], lty = 2, lwd = lwd)
        lines(times, Q.CI.upper[idx.lab[ii, 2], idx.lab[ii, 1],], lty = 2, lwd = lwd)
      }


      if(rug){
        pair.states = suStf$data[-nrow(suStf$data), '(state)'] == idx.lab[ii, 2] & suStf$data[-1, '(state)'] == idx.lab[ii, 1]
        same.person = suStf$data[-nrow(suStf$data), '(id)'] == suStf$data[-1, '(id)']
        act.tr.times = suStf$data[-1, suStf$tte][pair.states & same.person]
        rug(act.tr.times)
      }


    }
  }



  if(plot.Q.2d){

    if(length(smooth.2d.idx[[ii]]) == 0) warning('There are no 2D splines. plot.Q.2d = TRUE will be ignored.')

    if(is.na(match('view', names(match.call())))) view = object$msm.post.object$msm.post.object[[ii]]$smooth[[smooth.2d.idx[[ii]][1]]]$term # just choose the pair of the first 2d smooth if not provided
    if(is.na(match('theta', names(match.call())))) theta = 30
    if(is.na(match('phi', names(match.call())))) phi = 20
    if(is.na(match('cex.lab', names(match.call())))) cex.lab = 1.3

    times = newdata[, suStf$tte]
    which.plots = ifelse(is.null(which.plots), sum(suStf$whereQ != 0), which.plots)

    for(ii in 1:which.plots){

      if(length(smooth.2d.idx[[ii]]) > 0){

        if(is.null(cond.list.2d)) stop('The fixed variables for the plotting of the tensor interaction are missing. \nPlease provide them through the argument cond.list.2d.')

        obj = object$msm.post.object$msm.post.object[[ii]]

        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))

        vis.gam.msm(obj, view = view,
                    cond = cond.list.2d,
                    theta = theta, phi = phi, color = 'bw',
                    cex.lab = cex.lab, zlab = '',
                    s1 = idx.lab[ii, 2], s2 = idx.lab[ii, 1], ticktype = 'detailed',
                    type = 'link', zlim = c(0,0.8), xlab = 'Years since transplant', ...) # always link since we then add our exp function
        par(xpd = NA, srt = 98)  ## disable clipping and set string rotation
        text(-0.55, 0, paste(idx.lab[ii, 2], 'to', idx.lab[ii, 1], "transition intensity"), cex = 1.3)


      }



    }




  }






  list(Q.hist = Qmatr,
       Q.CI.lower.list = Q.CI.lower, Q.CI.upper.list = Q.CI.upper,
       full.X = full.X,
       # Qmatr.debug = Qmatr.debug, # for debugging/diagnostics
       Q.sim.hist = Q.CI.hist)


}















