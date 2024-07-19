# This function computes the trans probability matrix P given
# - a transition intensity matrix Q
# - the time lag
# - the way it's supposed to be computed (for now just eigendecomposition method)

# Returns:
# - P matrix
# - (perhaps also) the eigendecomposition elements, so can be used by derivative later on? Perhaps this is
# the solution to recycling the elements while keeping the computation of P separate from that of P deriv
# etc. ?
# Considerations on symmetry as in Jackson perhaps? Can this help for anything?


P.matr.comp = function(Qmatr, timelag, method = 'eigendecomp'){

  decomp.obj = NULL
  SS.obj = list(B.k.squaring.history = list(), B.k.scaling.history = list())


  if( method == 'eigendecomp'){

    tol = 15 #30 # 15 # tolerance for when we are checking whether the eigenvalues are distinct... too small, not small enough? Keep an eye on this

    decomp.obj = eigen(Qmatr)

    eigen.vec = decomp.obj$vectors
    eigen.val = decomp.obj$values

    # # **** Check for distinct eigenvalues ****
    # if( any(duplicated(round(eigen.val, tol))) ) stop('The eigenvalues of the Q matrix are not distinct, the eigendecompisition method cannot be used.')
    # # ****************************************

    eigen.vec.inv = solve(eigen.vec)

    # Save inverse in eigendecomp object, useful for later (dP part)
    decomp.obj$vectors.inv = eigen.vec.inv

    Pmatr = eigen.vec %*% diag( exp(eigen.val * timelag) )%*% eigen.vec.inv
  }


  if( method == 'analytic'){


    # allow user to define this tolerance? Does it create overflow issues? When?
    tol = 1e-15 # how small a difference counts as equality?

    Pmatr = matrix(0, nrow(Qmatr), ncol(Qmatr))

    # THIS HAS TO BE AUTOMATIZED !!!
    proc.type = 1 # i.e. an IDM
    ################################


    if(proc.type == 1){

      Pmatr[1,1] = exp(-(Qmatr[1,2] + Qmatr[1,3])*timelag)

      if( abs( (Qmatr[1,2] + Qmatr[1,3]) - Qmatr[2,3] ) > tol ){
        Pmatr[1,2] = Qmatr[1,2]/(Qmatr[1,2]+Qmatr[1,3]-Qmatr[2,3])*(exp(-Qmatr[2,3]*timelag) - exp(-(Qmatr[1,2] + Qmatr[1,3])*timelag))
      } else {
        Pmatr[1,2] = Qmatr[1,2]*timelag*exp(-(Qmatr[1,2] + Qmatr[1,3])*timelag)
      }

      Pmatr[1,3] = 1 - (Pmatr[1,1] + Pmatr[1,2])

      Pmatr[2,1] = 0
      Pmatr[2,2] = exp(-Qmatr[2,3]*timelag)
      Pmatr[2,3] = 1 - (Pmatr[2,1] + Pmatr[2,2])

      Pmatr[3,1] = 0
      Pmatr[3,2] = 0
      Pmatr[3,3] = 1

    } else if (proc.type == 2){ # not allowed yet since expressions for when intensities are equal are still missing


      # PROBABILITIES WITH STARTING STATE 1
      Pmatr[1,1] = exp( - Qmatr[1,2]*timelag )

      if(abs(Qmatr[2,3] -  Qmatr[1,2]) > tol){
        Pmatr[1,2] = -(exp(-Qmatr[2,3]*timelag) - exp(-Qmatr[1,2]*timelag))/(Qmatr[2,3] -  Qmatr[1,2])*Qmatr[1,2]
      } else {
        stop('This expression is not available yet (1-2).')
      }


      if(abs(Qmatr[3,4] - Qmatr[2,3]) > tol & abs(Qmatr[3,4] - Qmatr[1,2]) > tol & abs(Qmatr[2,3] - Qmatr[1,2]) > tol){
        Pmatr[1,3] = Qmatr[1,2] * Qmatr[2,3] / (Qmatr[2,3] - Qmatr[1,2]) *
          ( ( exp(-Qmatr[3,4]*timelag) - exp(-Qmatr[2,3]*timelag) )/(Qmatr[3,4] - Qmatr[2,3]) - ( exp(-Qmatr[3,4]*timelag) - exp(-Qmatr[1,2]*timelag) )/((Qmatr[3,4] - Qmatr[1,2])) )
      } else {
        stop('This expression is not available yet (1-3).')
      }


      if( abs(Qmatr[4,5] - Qmatr[3,4]) > tol & abs(Qmatr[4,5] - Qmatr[2,3]) > tol & abs(Qmatr[4,5] - Qmatr[1,2]) > tol & abs(Qmatr[3,4] - Qmatr[2,3]) > tol & abs(Qmatr[3,4] - Qmatr[1,2]) > tol & abs(Qmatr[2,3] - Qmatr[1,2]) > tol){
        Pmatr[1,4] = Qmatr[1,2] * Qmatr[2,3] * Qmatr[3,4] / (Qmatr[2,3] - Qmatr[1,2]) * (
          ( (exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[2,3]*timelag))/(Qmatr[4,5] - Qmatr[2,3]) - (exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[3,4]*timelag))/(Qmatr[4,5] - Qmatr[3,4])  )  / (Qmatr[3,4] - Qmatr[2,3])
            + ( (exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[3,4]*timelag))/(Qmatr[4,5] - Qmatr[3,4]) - (exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[1,2]*timelag))/(Qmatr[4,5] - Qmatr[1,2])  ) / (Qmatr[3,4] - Qmatr[1,2])
          )
      } else {
        stop('This expression is not available yet (1-4).')
      }


      if( abs(Qmatr[4,5] - Qmatr[3,4]) > tol & abs(Qmatr[4,5] - Qmatr[2,3]) > tol & abs(Qmatr[4,5] - Qmatr[1,2]) > tol & abs(Qmatr[3,4] - Qmatr[2,3]) > tol & abs(Qmatr[3,4] - Qmatr[1,2]) > tol & abs(Qmatr[2,3] - Qmatr[1,2]) > tol ){
        Pmatr[1,5] = Qmatr[1,2] * Qmatr[2,3] * Qmatr[3,4] * Qmatr[4,5] / (Qmatr[2,3] - Qmatr[1,2]) * (

          ( ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[3,4]*timelag) - 1)/Qmatr[3,4] )/(Qmatr[4,5] - Qmatr[3,4])
            - ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[2,3]*timelag) - 1)/Qmatr[2,3] )/(Qmatr[4,5] - Qmatr[2,3])  ) / (Qmatr[3,4] - Qmatr[2,3])

          + ( ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[1,2]*timelag) - 1)/Qmatr[1,2] )/(Qmatr[4,5] - Qmatr[1,2])
              - ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[3,4]*timelag) - 1)/Qmatr[3,4] )/(Qmatr[4,5] - Qmatr[3,4])  ) / (Qmatr[3,4] - Qmatr[1,2])
          )
      } else {
        stop('This expression is not available yet (1-5).')
      }



      # PROBABILITIES WITH STARTING STATE 2
      Pmatr[2,2] = exp(-Qmatr[2,3]*timelag)

      if(abs(Qmatr[3,4] -  Qmatr[2,3]) > tol){
        Pmatr[2,3] = -(exp(-Qmatr[3,4]*timelag) - exp(-Qmatr[2,3]*timelag))/(Qmatr[3,4] -  Qmatr[2,3])*Qmatr[2,3]
      } else {
        stop('This expression is not available yet (2-3).')
      }

      if(abs(Qmatr[4,5] - Qmatr[3,4]) > tol & abs(Qmatr[4,5] - Qmatr[2,3]) > tol & abs(Qmatr[3,4] - Qmatr[2,3]) > tol){
        Pmatr[2,4] = Qmatr[2,3] * Qmatr[3,4] / (Qmatr[3,4] - Qmatr[2,3]) *
          ( ( exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[3,4]*timelag) )/(Qmatr[4,5] - Qmatr[3,4]) - ( exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[2,3]*timelag) )/((Qmatr[4,5] - Qmatr[2,3])) )
      } else {
        stop('This expression is not available yet (2-4).')
      }


      if(abs(Qmatr[4,5] - Qmatr[3,4]) > tol & abs(Qmatr[4,5] - Qmatr[2,3]) > tol & abs(Qmatr[3,4] - Qmatr[2,3]) > tol){
        Pmatr[2,5] = Qmatr[2,3] * Qmatr[3,4] * Qmatr[4,5] / (Qmatr[3,4] - Qmatr[2,3]) * (
          ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[2,3]*timelag) - 1)/Qmatr[2,3] )/(Qmatr[4,5] - Qmatr[2,3])
          - ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[3,4]*timelag) - 1)/Qmatr[3,4] )/(Qmatr[4,5] - Qmatr[3,4])
        )
      } else {
        stop('This expression is not available yet (2-5).')
      }


      # PROBABILITIES WITH STARTING STATE 3
      Pmatr[3,3] = exp(-Qmatr[3,4]*timelag)

      if(abs(Qmatr[4,5] - Qmatr[3,4]) > tol){
        Pmatr[3,4] = - Qmatr[3,4] * (exp(-Qmatr[4,5]*timelag) - exp(-Qmatr[3,4]*timelag) ) / (Qmatr[4,5] - Qmatr[3,4])
      } else {
       stop('This expression is not available yet (3-4).')
      }

      if(abs(Qmatr[4,5] - Qmatr[3,4]) > tol){
      Pmatr[3,5] = Qmatr[3,4] * Qmatr[4,5] / (Qmatr[4,5] - Qmatr[3,4]) * ( (exp(-Qmatr[4,5]*timelag) - 1)/Qmatr[4,5] - (exp(-Qmatr[3,4]*timelag) - 1)/Qmatr[3,4] )
      } else {
        stop('This expression is not available yet (3-5).')
      }


      # PROBABILITIES WITH STARTING STATE 4
      Pmatr[4,4] = exp(-Qmatr[4,5]*timelag)
      Pmatr[4,5] = 1 - exp(-Qmatr[4,5]*timelag)

      # PROBABILITIES WITH STARTING STATE 5
      Pmatr[5,5] = 1

    }





  }




  if( method == 'scaling&squaring'){

    warning('The scaling&squaring method is still a work in progress. Please use the other methods instead as they are faster.')

    # Values of N and m proposed by Zhong & Williams (1994)
    N = 20
    m = 4

    # Matrix setup
    A = Qmatr*timelag
    Id.m = diag(dim(A)[1])


    if( 1/2^N * max(abs(A)) >= 1  ){
      N = ceiling(log2(max(abs(A))))
      write(paste('Largest Q value found was', max(abs(A)), 'so N was set to', N, '\n'), file = 'testing_adaptive_S&S.txt', append = TRUE)
    }

    # *** Computation of B^[N] *** as proposed in Young (2004), with correction of typos
    # (check handwritten notes for this)
    # NOTE: we only work on the terms of the taylor expansion from the
    # (k = 1)th and on, "to avoid the unnecessary errors otherwise
    # caused by T_a being very small" (Zhong & Williams, 1994).

    B.k.min.1 = 1/2^N * A # summation starts at k = 1 since we skip the first term of the expansion, i.e. the identity matrix
    B.N.tay.expansion = B.k.min.1

    # We set this up to recycle quantities already computed for the computation of first and second derivatives
    SS.obj$B.k.squaring.history[[1]] = B.k.min.1

    for(k in 2:m){ # start addition of following terms from k = 2, since k = 1 is starting point
      B.k = 1/k * 1/2^N * B.k.min.1 %*% A
      B.N.tay.expansion = B.N.tay.expansion + B.k
      B.k.min.1 = B.k

      # Save this to recycle quantities
      SS.obj$B.k.squaring.history[[k]] = B.k.min.1

    }

    # # Save also B^[N] always with the purpose of recycling quantities ---- ACTUALLY WE DO NOT NEED THIS since we have full history
    # SS.obj$B.N = B.N.tay.expansion

    # *** Computation of B^[0] *** (according to recursive formula proposed in Young (2004) and Zhong & Williams (1994))
    B.k = B.N.tay.expansion
    # Save this to recycle quantities
    SS.obj$B.k.scaling.history[[1]] = B.k

    for(k in N:1){
      B.k = 2 * B.k + B.k %*% B.k
      # Save this to recycle quantities
      SS.obj$B.k.scaling.history[[N-k+2]] = B.k # +2 because these quantities go from N to 0 (so in code we need to shift to N+1 to 1)
    }

    # # Finally, save B^[0] always with the purpose of recycling quantities ---- ACTUALLY WE DO NOT NEED THIS since we have full history
    # SS.obj$B.0 = B.k

    # Finally, the matrix exponential follows from formula: B^[0] + I = exp(A)
    Pmatr = B.k + Id.m


  }

  list(Pmatr = Pmatr,
       decomp.obj = decomp.obj,
       SS.obj = SS.obj)

}
