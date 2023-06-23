
#' @title Predict and plot the transition probabilities
#'
#' @description
#' Function to predict and plot the estimated transition probabilities (and confidence intervals).
#'
#' @param object Fitted model object.
#' @param newdata Dataframe containing the profile for which one wished to obtain the predicted transition probabilities.
#' @param get.CI Whether to compute the confidence intervals.
#' @param n.sim.CI Number of simulations to be used for confidence intervals computation.
#' @param prob.lev Probability level of confidence intervals.
#' @param plot.P Whether to output plots of transition probabilities.
#' @param which.plots Number between 1 and the maximum number of non-null transition probabilities. This can be used if only some plots are to be plotted.
#' @param rug Whether to include a rugplot of the observed transition times.
#' @param ... Other graphical arguments.
#'
#' @usage P.pred(object, newdata, get.CI = TRUE,
#'        n.sim.CI = 1000, prob.lev = 0.05,
#'        plot.P = FALSE, which.plots = NULL,
#'        rug = FALSE, ...)
#'
#' @return Estimated transition probabilities (and confidence intervals).
#' \item{P.pred}{Predicted transition probability matrix corresponding to the time horizon specified in \code{newdata}. This is a \code{nstates x nstates} matrix.}
#' \item{P.CI.lower}{Matrix containing the lower bounds of the confidence intervals for the predicted transition probability matrix.}
#' \item{P.CI.upper}{Matrix containing the upper bounds of the confidence intervals for the predicted transition probability matrix.}
#' \item{P.hist}{List of predicted transition probability matrices computed at each time point specified in \code{newdata}. This is a \code{nstates x nstates x n.pred} array, where \code{n.pred} is the number of rows in \code{newdata}.}
#' \item{P.CI.lower.hist}{List of matrices containing the lower bounds of the confidence intervals for each predicted transition probability matrix in \code{P.hist}..}
#' \item{P.CI.upper.hist}{List of matrices containing the upper bounds of the confidence intervals for each predicted transition probability matrix in \code{P.hist}.}
#' \item{full.X}{Full design matrix corresponding to the \code{newdata} provided.}
#' \item{P.sim.hist}{List of transition probability matrices simulated to obtain the confidence intervals at each time point from \code{newdata}. May be useful to quickly obtain intervals for a different confidence level.}
#'
#'
#' @export
#'
#' @seealso \code{\link{fmsm}}
#'

P.pred = function(object, newdata,
                  get.CI = TRUE, n.sim.CI = 1000,
                  prob.lev = 0.05, plot.P = FALSE,
                  # plot.P.stacked = FALSE,
                  which.plots = NULL, rug = FALSE, ...){


  # if(plot.P & plot.P.stacked) stop('For plotting, you have to choose one between the traditional and stacked visualisation.')
  # if(plot.P.stacked & get.CI) warning('Confidence intervals will not be plotted with stacked visualisation.')

  P.CI.hist = P.CI.lower = P.CI.upper = discard = NULL

  params.hat = object$msm.fit.object$fit$argument
  sp.hat = object$msm.fit.object$sp

  n.pred = nrow(newdata)
  if(n.pred < 2) stop('newdata needs to have at least two rows.')

  suStf = object$suStf
  pmethod = suStf$pmethod
  death = suStf$death

  full.X = c()
  for(i in 1:length(object$msm.post.object$msm.post.object)){
    mod = object$msm.post.object$msm.post.object[[i]]
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

  timelags = diff(newdata[,suStf$tte])

  P.hist = array(dim = c(suStf$nstates, suStf$nstates, n.pred))
  P.hist[,,1] = P.pred = diag(rep(1, suStf$nstates)) # starting point is the identity matrix (this is P(0, 0) = I)

  for(i in 1:(n.pred-1)){
    P.pred = P.pred %*% P.matr.comp(Qmatr[,,i], timelags[i])$Pmatr
    P.hist[,,i+1] = P.pred
  }

  idx.lab = which(t(P.pred) != 0, arr.ind = T)


  if(get.CI){

    # Confidence intervals ####
    beta.mu = params.hat
    beta.sigma = matrix(object$msm.post.object$Vb, nrow = length(beta.mu), ncol = length(beta.mu))
    # set.seed(24)
    bs = GJRM::rMVN(n = n.sim.CI, mean = beta.mu, sigma = beta.sigma)

    P.CI.hist = array(dim = c(suStf$nstates, suStf$nstates, n.pred, n.sim.CI))
    discard = c()

    for(i.CI in 1:n.sim.CI){

      Qmatr = Q.matr.setup.general(params = bs[i.CI,][suStf$pos.optparams], nstates = suStf$nstates,
                                   full.X = full.X, start.pos.par = suStf$start.pos.par,
                                   l.short.formula = suStf$l.short.formula, whereQ = suStf$whereQ,
                                   firstD = FALSE, secondD = FALSE, bound.eta = FALSE,
                                   pos.optparams = suStf$pos.optparams, pos.optparams2 = suStf$pos.optparams2)$Qmatr

      timelags = diff(newdata[,suStf$tte])

      P.CI.hist[, , 1, i.CI] = P.pred = diag(rep(1, suStf$nstates))

      for(i in 1:(n.pred-1)){
        P.temp = try(P.matr.comp(Qmatr[,,i], timelags[i])$Pmatr)

        if(!is.null(attr(P.temp, 'class'))){
          discard = c(discard, i.CI)
          stop('Something is going wrong in the computation of the P matrix confidence intervals. Contact author for help.')
          break
        }

        P.pred = P.pred %*% P.matr.comp(Qmatr[,,i], timelags[i])$Pmatr
        P.CI.hist[, , i+1, i.CI] = P.pred
      }
      rm(P.temp)

    }


    P.CI.lower = P.CI.upper = array(0, dim = c(suStf$nstates, suStf$nstates, n.pred))

    for(ii in 1:sum(P.pred != 0)){
      P.CI.lower[idx.lab[ii, 2], idx.lab[ii, 1],] = matrixStats::rowQuantiles(P.CI.hist[idx.lab[ii, 2], idx.lab[ii, 1],,], probs = c(prob.lev/2), na.rm = TRUE)
      P.CI.upper[idx.lab[ii, 2], idx.lab[ii, 1],] = matrixStats::rowQuantiles(P.CI.hist[idx.lab[ii, 2], idx.lab[ii, 1],,], probs = c(1 - prob.lev/2), na.rm = TRUE)
    }

    # if(i.CI/n.sim.CI*100 %% 20 == 0) print(paste(i.CI/n.sim.CI*100, '% of simulations for the CIs completed.'))

  }


  if(plot.P){

    if(is.na(match('xlab', names(match.call())))) xlab = suStf$tte
    if(is.na(match('cex.lab', names(match.call())))) cex.lab = 1.2
    if(is.na(match('cex.axis', names(match.call())))) cex.axis = 1.25

    lwd.P = 1
    lwd.CI = 1

    times = newdata[, suStf$tte]
    which.plots = ifelse(is.null(which.plots), sum(P.pred != 0), which.plots)

    for(ii in 1:which.plots){

      if(death & idx.lab[ii, 2] == suStf$nstates & idx.lab[ii, 1] == suStf$nstates) break

      Ps = P.hist[idx.lab[ii, 2], idx.lab[ii, 1],]

      plot(times, Ps, type = 'l', lwd = lwd.P,
           ylim = c(0,1),
           xlab = '',
           ylab = paste(idx.lab[ii, 2], 'to', idx.lab[ii, 1], 'trans. probability'),
           cex.lab = cex.lab, cex.axis = cex.axis)
      title(xlab = xlab, line = 4, cex.lab = cex.lab)


      if(get.CI) {
        lines(times, P.CI.lower[idx.lab[ii, 2], idx.lab[ii, 1],], lty = 2, lwd = lwd.CI)
        lines(times, P.CI.upper[idx.lab[ii, 2], idx.lab[ii, 1],], lty = 2, lwd = lwd.CI)
      }

      if(rug){
        pair.states = suStf$data[-nrow(suStf$data), '(state)'] == idx.lab[ii, 2] & suStf$data[-1, '(state)'] == idx.lab[ii, 1]
        same.person = suStf$data[-nrow(suStf$data), '(id)'] == suStf$data[-1, '(id)']
        act.tr.times = suStf$data[-1, suStf$tte][pair.states & same.person]
        rug(act.tr.times)
      }


    }



   }
# else if (plot.P.stacked){
#
#     times = newdata[, suStf$tte]
#
#     for(ii in 1:suStf$nstates){
#
# # SOMETHING - NEED TO COMPLETE THIS
#
#     }
#
#     for(ii in 1:sum(P.pred != 0)){
#
#       Ps = P.hist[idx.lab[ii, 2], idx.lab[ii, 1],]
#
#       plot(times, Ps, type = 'l', lwd = lwd.P,
#            ylim = c(0,1),
#            # xlab = '',
#            ylab = paste(idx.lab[ii, 2], 'to', idx.lab[ii, 1], 'trans. probability'),
#            # xaxt = 'n',
#            cex.lab = 1.4, cex.axis = 1.25)
#
#       plot(times, probs.test$pstate1, type = 'l', lwd = 2, col = 'grey50', ylim = c(0,1),
#            xlab = 'Follow-up time', main = main, ylab = ylab)
#       polygon(c(times, rev(times)), c(probs.test$pstate1+probs.test$pstate2, probs.test$pstate1+probs.test$pstate2+probs.test$pstate3),
#               col = 'grey90',
#               border = NA)
#       lines(times, probs.test$pstate1+probs.test$pstate2, lwd = 2, col = 'grey70')
#       polygon(c(times, rev(times)), c(probs.test$pstate1, rev(probs.test$pstate1+probs.test$pstate2)),
#               col = 'grey70',
#               border = NA)
#       polygon(c(rev(times), times), c(rep(0, length(times)), probs.test$pstate1),
#               col = 'grey50',
#               border = NA)
#       lines(times, probs.test$pstate1+probs.test$pstate2+probs.test$pstate3, lwd = 2, col = 'grey90')
#
#
#
#     }
#
#
#
#   }







  list(P.pred = P.hist[,,dim(P.hist)[3]],
       P.CI.lower = P.CI.lower[,,n.pred], P.CI.upper = P.CI.upper[,,n.pred],
       P.hist = P.hist,
       P.CI.lower.hist = P.CI.lower, P.CI.upper.hist = P.CI.upper,
       # discard = discard,
       full.X = full.X,
       # Qmatr.debug = Qmatr.debug, # for debugging/diagnostics
       # Qmatr = Qmatr,
       P.sim.hist = P.CI.hist)


}















