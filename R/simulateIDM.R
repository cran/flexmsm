
#' Function to predict and plot the estimated transition intensities (and confidence intervals).
#'
#' @param N Total number of individuals.
#' @param seed Seed used for the simulation.
#' @param og.12 If \code{TRUE} a common shape for the first transition is used. If \code{FALSE} a sinusoid is used.
#'
#' @return Simulated data generated from an Illness-Death model (IDM).
#'
#' @export
#'

simulateIDM = function(N = N, seed = seed, og.12 = TRUE){

  set.seed(seed = seed)

  # # **** NO COVARIATES FOR NOW ****
  # # Set linear covariate effects values
  # beta.dage = 0.015
  # beta.pdiag = 0.2
  #
  # covs = cbind(dage.sim, IHD.sim)
  # beta.lin = rbind(beta.dage, beta.pdiag)

  exp.covs = 1 #exp(covs %*% beta.lin)

  if(og.12){
    # T.12
    mu = 1.25
    sigma = 1

    u  <- runif(N)
    val <- 1-exp( log( 1-u )*exp.covs ) # more complex version needed rather than just u because of covariates
    sim.T.12 <- qlnorm(val , mu, sigma)
  } else { # this is simulating the 1->2 time-to-event using the cosine shape
    # q.12 = function(t) 0.1 * cos((2*pi/10) * t + pi/2) + 0.15
    H.12 = function(t) 0.15 * t + 1/(2*pi) * sin((2*pi/10) * t + pi/2) - 1/(2*pi)

    u <- runif(N)
    vals <- -log(1-u)
    maxfup <- 50 # maximum follow-up time

    pit.fun = function(t, val) H.12(t) - val # Probability Integral Transform function

    ind.cens <- (vals > H.12(maxfup))
    sim.T.12 <- rep(0, N) # simulated times
    sim.T.12[ind.cens] <- maxfup # censored times
    ind.obs <- (1:N)[!ind.cens] # observed times
    # Probability Integral Transform based on the cumulative hazard
    for(i in ind.obs) sim.T.12[i] <- uniroot(f = pit.fun, interval = c(0, maxfup), val = vals[i], maxiter = 10000)$root
  }


  # T.13
  lambda = exp(-2.5)

  u  <- runif(N)
  val <- 1-exp( log( 1-u )*exp.covs ) # more complex version needed rather than just u because of covariates
  sim.T.13 <- qexp(val, rate = lambda)


  # T.23 (for this we need conditional distribution - Ardo book Sec. 2.3.3)
  rate = exp(-2.5)
  shape = 0.1

  u  <- runif(N)
  val <- 1-exp( log( 1-u ) ) # just u because in conditional case time-fixed covariates disappear
  T.gomp.cond = function(x, t, a = shape, b = rate) (log(b*exp(a*t) - a*log(x)) - log(b))/a
  sim.T.23 <- T.gomp.cond(x = val, t = sim.T.12, a = shape, b = rate)


  # impose censoring
  t.obs.times = 0:15

  time = id = state = c()

  for(i in 1:N){

    if(min(sim.T.12[i], sim.T.13[i]) == sim.T.12[i]){

      upper.t = ceiling(sim.T.12[i])
      upper.t.2 = ceiling(sim.T.23[i])

      if(sim.T.12[i] > 15){
        time = c(time, t.obs.times)
        state = c(state, rep(1, length(t.obs.times)))
        id = c(id, rep(i, length(t.obs.times)))

        next
      }

      if(sim.T.23[i] > 15){
        time = c(time, t.obs.times)
        state = c(state, rep(1, length(t.obs.times[1:upper.t])),
                  rep(2, length(t.obs.times[(upper.t+1):length(t.obs.times)])))
        id = c(id, rep(i, length(t.obs.times)))

        next
      }


      if(upper.t == upper.t.2){ # trans 1->2 + 2->3 happens between observations times so only 1->3 is "actually observed"
        time = c(time, t.obs.times[1:upper.t.2], sim.T.23[i])
        state = c(state, rep(1, length(t.obs.times[1:upper.t.2])), 3)
        id = c(id, rep(i, length(t.obs.times[1:upper.t.2])+1))
      } else {
        time = c(time, t.obs.times[1:upper.t.2], sim.T.23[i])
        state = c(state, rep(1, length(t.obs.times[1:upper.t])),
                  rep(2, length(t.obs.times[(upper.t+1):upper.t.2])), 3)
        id = c(id, rep(i, length(t.obs.times[1:upper.t.2])+1))
      }

    } else {

      if(sim.T.13[i] > 15){
        time = c(time, t.obs.times)
        state = c(state, rep(1, length(t.obs.times)))
        id = c(id, rep(i, length(t.obs.times)))

        next
      }

      upper.t = ceiling(sim.T.13[i])
      time = c(time, t.obs.times[1:upper.t], sim.T.13[i])
      state = c(state, rep(1, length(t.obs.times[1:upper.t])), 3)
      id = c(id, rep(i, length(t.obs.times[1:upper.t])+1))
    }

    # if( length(id) != length(state) | length(id) != length(time) | length(state) != length(time)) print(i)

  }

  dataSim = data.frame(cbind(id, time, state))
  names(dataSim) = c('PTNUM', 'years', 'state')



  dataSim
}
