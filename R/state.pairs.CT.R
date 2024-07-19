

#' Function to extract state pair counts and observed (right) times.
#'
#' @param formula Model specification.
#' @param data Data.
#' @param whereQ Placement of allowed transition intensities. Only for internal use. Defaults to \code{NULL} and is obtained automatically when \code{formula} is provided.
#' @param nstates Total number of states. Only for internal use. Defaults to \code{NULL} and is obtained automatically when \code{formula} is provided.
#' @param time Name of variable containing the time-to-event.
#' @param state Name of variable containing the states.
#' @param id Name of variable containing the unique code identifying the individuals.
#'
#' @return A table with the state-pair counts and a list with the observed (right) times for each transition.
#'
#' @export
#'

state.pairs.CT = function(formula = NULL, data = NULL, whereQ = NULL, nstates = NULL,
                          time = NULL, state = NULL, id = NULL){ # pass these on as strings

  if((!is.null(formula) & !is.null(whereQ)) | (is.null(formula) & is.null(whereQ))) stop('You need to provide either formula or whereQ and nstates.')


  if(!is.null(formula)){

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
  }


  counts = matrix(0, nrow = nstates, ncol = nstates)
  act.tr.times = list()
  act.who = list()
  idx.lab = which(t(whereQ) != 0, arr.ind = T)

  for(i in 1:sum(whereQ != 0)){
    pair.states = data[-nrow(data), state] == idx.lab[i, 2] & data[-1, state] == idx.lab[i, 1]
    same.person = data[-nrow(data), id] == data[-1, id]
    act.who[[i]] = which(pair.states & same.person)
    act.tr.times[[i]] = data[-1, time][pair.states & same.person]

    counts[idx.lab[i, 2], idx.lab[i, 1]] = length(act.tr.times[[i]])
  }

  # And now the diagonal elements
  for(i in 1:nstates){
    pair.states = data[-nrow(data), state] == i & data[-1, state] == i
    same.person = data[-nrow(data), id] == data[-1, id]
    act.who[[length(act.who)+1]] = which(pair.states & same.person)
    act.tr.times[[length(act.tr.times)+1]] = data[-1, time][pair.states & same.person]

    counts[i, i] = length(act.tr.times[[length(act.tr.times)]])
  }


  rownames(counts) = colnames(counts) = paste('state', 1:nstates)

  tmp.pairs = cbind(data[-nrow(data), state], data[-1, state])[same.person,]
  full.counts = table(tmp.pairs[,1], tmp.pairs[,2])

  list(counts = counts,
       act.tr.times = act.tr.times,
       act.who = act.who,
       full.counts = full.counts)

}



