#' Internal function
#'
#' Internal function needed for setup of Q matrix and its first and second derivative.
#'
#' @param params Parameters vector.
#' @param nstates Number of states.
#' @param full.X Full design matrix.
#' @param start.pos.par Positions within full parameters vector of starting point for each sub-parameters vector corresponding to each transition intensity specification.
#' @param l.short.formula Number of transitions.
#' @param whereQ Positions within Q matrix of not-null transition intensities.
#' @param firstD Whether the first derivative of the Q matrix should be computed.
#' @param secondD Whether the second derivative of the Q matrix should be computed.
#' @param bound.eta Whether to bound the additive predictor, defaults to \code{FALSE}. This is only used for debugging purposes, do not change.
#' @param pos.optparams Vector with positions of parameters vector in the form used by the optimization algorithm (i.e. when one or more parameters are constrained to be equal these will only appear once).
#' @param pos.optparams2 Like \code{pos.optparams} but the count is not stopped at the constrained parameters.
#'
#' @return Q matrix and its first and second derivatives with respect to the parameters vector.
#' @export
#'
Q.matr.setup.general = function(params, nstates, full.X, start.pos.par, l.short.formula, whereQ, firstD = TRUE, secondD = TRUE, bound.eta = FALSE,
                                    pos.optparams, pos.optparams2){ # last argument added as of 05/10/2022 since we are generalising code



  comp.par.mapping = NULL

  # If could have unique global formula I could avoid repeating the dataframes and have something like Jackson (with
  # corresponding parameters set to 0 if not present for given transition)... perhaps this approach is neater ?
  for(i in 1:l.short.formula){

    matr = full.X[,start.pos.par[i]:(start.pos.par[i+1]-1)]

    if(bound.eta){
      eta =  matr %*% params[start.pos.par[i]:(start.pos.par[i+1]-1)]
      eta[eta > 1] = 1
    } else {
      eta =  matr %*% params[start.pos.par[i]:(start.pos.par[i+1]-1)]
    }

    # q matrix values
    assign(paste('Qvals', i, sep = ''), exp( eta ) )

  }


  # !!!! ************************************ !!!! #
  # CHECK CORRECTNESS OF QMATR, DQMATR AND D2QMATR #
  # !!!! ************************************ !!!! #

  Qmatr = array(0, dim = c(nstates, nstates, nrow(full.X)))
  for(i in 1:nstates){
    for(j in 1:nstates){
      if(whereQ[i,j] != 0) Qmatr[i, j, ] = eval(as.name(paste('Qvals', whereQ[i,j], sep = '')))
    }
    # Insert diagonal elements once all of the elements of a given row have been taken care of
    Qmatr[i, i, ] = - apply(matrix(Qmatr[i, , ], ncol = dim(Qmatr)[3]), 2, sum)
  }







  if(firstD){

    dQmatr = array(0, dim = c(nstates, nstates, length(unique(pos.optparams)), nrow(full.X))) # note: on 17/02/2022 changed params to optparams here to allow and account for constraining of parameters

    for(i in 1:nstates){
      for(j in 1:nstates){
        if(whereQ[i,j] != 0){
          for(l in start.pos.par[whereQ[i,j]]:(start.pos.par[whereQ[i,j]+1]-1)) dQmatr[i, j, pos.optparams[l], ] = Qmatr[i, j, ] * full.X[,l] # on 17/02/2022 changed indexing in dQmatr when assigning to dQmatr to account for constrained parameters
        }
      }
      # Insert diagonal elements once all of the elements of a given row have been taken care of
      dQmatr[i, i, , ] = - apply(array(dQmatr[i, , , ], dim = c(dim(dQmatr)[2:4])), c(2,3), sum)
    }





  } else {

    dQmatr = NULL

  }






  if(secondD){

    start.pos.optparams = c(pos.optparams[start.pos.par[-length(start.pos.par)]], max(pos.optparams)+1) # remember always need index right after the last as well (so everything is correct when taking differences of positions to compute the number of parameters)
    actual.n.2param = sum(diff(start.pos.optparams)*(diff(start.pos.optparams)+1)*0.5)

    # If there are constrained parameters we need to add some derivatives (the first piece only
    # accounts for derivative within the parameters of one specific transition, removing the constrained parameters)
    for(idx in 1:(length(start.pos.par)-1)){
      constr.idx = pos.optparams2[start.pos.par[idx]:(start.pos.par[idx+1]-1)]
      unconstr.idx = start.pos.par[idx]:(start.pos.par[idx+1]-1)
      if(any(unconstr.idx != constr.idx)){
        actual.n.2param = actual.n.2param + sum(unconstr.idx != constr.idx) * sum(unconstr.idx == constr.idx)
      }
    }

    d2Qmatr = array(0, dim = c(nstates, nstates, actual.n.2param, nrow(full.X)))
    comp.par = 1 # compact parameter counter

    comp.par.mapping = matrix(NA, nrow = length(unique(pos.optparams)), ncol = length(unique(pos.optparams)))

    for(i in 1:nstates){
      for(j in 1:nstates){
        if(whereQ[i,j] != 0){
          for(l in start.pos.par[whereQ[i,j]]:(start.pos.par[whereQ[i,j]+1]-1)){
            for(k in l:(start.pos.par[whereQ[i,j]+1]-1)) {

              if( is.na(comp.par.mapping[min(pos.optparams[l], pos.optparams[k]), max(pos.optparams[l], pos.optparams[k])]) ){
                d2Qmatr[i, j, comp.par, ] = dQmatr[i, j, pos.optparams[k], ] * full.X[,l]
                comp.par.mapping[min(pos.optparams[l], pos.optparams[k]), max(pos.optparams[l], pos.optparams[k])] = comp.par
                comp.par = comp.par + 1
              } else {
                d2Qmatr[i, j, comp.par.mapping[min(pos.optparams[l], pos.optparams[k]), max(pos.optparams[l], pos.optparams[k])], ] = dQmatr[i, j, pos.optparams[k], ] * full.X[,l]
              }


            }
          }
        }
      }
      # Insert diagonal elements once all of the elements of a given row have been taken care of
      d2Qmatr[i, i, , ] = - apply(array(d2Qmatr[i, , , ], dim = c(dim(d2Qmatr)[2:4])), c(2,3), sum)
      endd2 = Sys.time()
    }
    end = Sys.time()


    if( (comp.par - 1) != actual.n.2param ) warning(paste("Something is not right in d2Q mapping come and take a look! comp.par is", comp.par, "while actual.n.2param is", actual.n.2param))

    # # *************************************


  } else {

    d2Qmatr = NULL

  }



  list(Qmatr = Qmatr, dQmatr = dQmatr, d2Qmatr = d2Qmatr,
       params = params,
       comp.par.mapping = comp.par.mapping)

}
