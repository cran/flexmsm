# # This function computes the derivative of trans probability matrix P given
# - a transition intensity matrix Q
# - the time lag
# - the derivative of the Q matrix

# NOTE: decomp.obj.i is needed for eigendecomp method while Pmatr. and Qmatr.i are needed for analytical method.

# For now implemented only for method = 'eigendecomp'.
#
#
dP.matr.comp = function(nstates, l.params, timelag, dQmatr.i, method = 'eigendecomp', decomp.obj.i = NULL, SS.obj.i = NULL, Pmatr.i = NULL, Qmatr.i = NULL){

  if(method == 'eigendecomp'){

    tol = 30 #15 # tolerance for when we are checking whether the eigenvalues are distinct... too small, not small enough? Keep an eye on this

    eigen.vec.i = decomp.obj.i$vectors
    eigen.val.i = decomp.obj.i$values

    # **** Check for distinct eigenvalues **** perhaps not needed since should get blocked at the level of P.matr.comp() function (which is always called before)
    if( any(duplicated(round(eigen.val.i, tol))) ) stop('The eigenvalues of the Q matrix are not distinct, the eigendecompisition method cannot be used.')
    # ****************************************

    eigen.vec.inv.i = decomp.obj.i$vectors.inv

    # Initialize matrix
    dPmatr.i = array(0, dim = c(nstates, nstates, l.params)) # because it's for single individual


    # *** NEW! Setup of required quantities (outside of for loop) ***************
    eem1 = matrix(exp(rep(eigen.val.i, nstates)*timelag), ncol = nstates)
    eem2 = t(eem1)
    diag(eem2) = 0

    num = eem1 - eem2

    eem1 = matrix(rep(eigen.val.i, nstates), ncol = nstates)
    diag(eem1) = 1/timelag
    eem2 = t(eem1)
    diag(eem2) = 0

    denom = eem1 - eem2
    # ****************************************************************************

    for(r in 1:l.params){

      G.r = eigen.vec.inv.i %*% dQmatr.i[,,r] %*% eigen.vec.i

      # # Old code (inside the for loop)
      # eem1 = matrix(exp(rep(eigen.val.i, nstates)*timelag), ncol = nstates)
      # eem2 = t(eem1)
      # diag(eem2) = 0
      #
      # num = eem1 - eem2
      #
      # eem1 = matrix(rep(eigen.val.i, nstates), ncol = nstates)
      # diag(eem1) = 1/timelag
      # eem2 = t(eem1)
      # diag(eem2) = 0
      #
      # denom = eem1 - eem2

      dPmatr.i[,,r] = eigen.vec.i %*% (G.r * num / denom) %*% eigen.vec.inv.i

    }


    # ******
    SS.obj.i = NULL # NOT SURE IF THIS IS NEEDED, WILL IT JUST TAKE NULL AUTOMATICALLY?

  }



  if( method == 'analytic'){

    # allow user to define this tolerance? Does it create overflow issues? When?
    tol = 1e-15 # how small a difference counts as equality?

    dPmatr.i = array(0, dim = c(nstates, nstates, l.params))

    # Repeated quantities ***
    rq1 = exp(-Qmatr.i[2,3]*timelag) - exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)
    drq1 = - timelag * dQmatr.i[2,3,] * exp(-Qmatr.i[2,3]*timelag) + timelag * (dQmatr.i[1,2,] + dQmatr.i[1,3,]) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)

    rq2 = Qmatr.i[1,2] + Qmatr.i[1,3] - Qmatr.i[2,3]
    drq2 = dQmatr.i[1,2,] + dQmatr.i[1,3,] - dQmatr.i[2,3,]
    # *****

    dPmatr.i[1,1,] = Pmatr.i[1,1] * timelag * dQmatr.i[1,1,]

    if( abs( (Qmatr.i[1,2] + Qmatr.i[1,3]) - Qmatr.i[2,3] ) > tol ){
      dPmatr.i[1,2,] = dQmatr.i[1,2,] / rq2 * rq1 - Qmatr.i[1,2] * drq2 / rq2^2 * rq1 + Qmatr.i[1,2] / rq2 * drq1
    } else {
      dPmatr.i[1,2,] = exp(-(Qmatr.i[1,2]+Qmatr.i[1,3])*timelag) * timelag * (dQmatr.i[1,2,] - timelag * Qmatr.i[1,2]*(dQmatr.i[1,2,]+dQmatr.i[1,3,]))
    }

    dPmatr.i[1,3,] = - (dPmatr.i[1,1,] + dPmatr.i[1,2,])

    dPmatr.i[2,1,] = 0
    dPmatr.i[2,2,] = -timelag * Pmatr.i[2,2] * dQmatr.i[2,3,]
    dPmatr.i[2,3,] = - (dPmatr.i[2,1,] + dPmatr.i[2,2,])

    dPmatr.i[3,1,] = 0
    dPmatr.i[3,2,] = 0
    dPmatr.i[3,3,] = 0

    # ******
    SS.obj.i = NULL # NOT SURE IF THIS IS NEEDED, WILL IT JUST TAKE NULL AUTOMATICALLY?

  }




  if( method == 'scaling&squaring'){

    warning('The scaling&squaring method is still a work in progress. Please use the other methods instead as they are faster.')


    # Values of N and m proposed by Zhong & Williams (1994)
    N = 20
    m = 4

    # Matrix setup
    A = Qmatr.i*timelag
    dA = dQmatr.i*timelag
    Id.m = diag(dim(A)[1])

    if( 1/2^N * max(abs(A)) >= 1  ){
      N = ceiling(log2(max(abs(A))))
      write(paste('Largest Q value found was', max(abs(A)), 'so N was set to', N, '\n'), file = 'testing_adaptive_S&S.txt', append = TRUE)
    }


    # Initialize matrix
    dPmatr.i = array(0, dim = c(nstates, nstates, l.params)) # because it's for single individual

    # Expand SS.obj to recycle first derivative quantities for the second derivatives computation
    SS.obj.i$d.histories = list()

    for(r in 1:l.params){

      # Initializes the rth history list
      SS.obj.i$d.histories[[r]] = list(dB.k.squaring.history = list(), dB.k.scaling.history = list())


      # *** Computation of B_{,r}^[N] *** as proposed in Young (2004), with correction of typos
      # B.k.min.1 = 1/2^N * A
      dB.k.min.1.r = 1/2^N * dA[,,r] # k = 1 term of summation

      dB.N.tay.expansion.r = dB.k.min.1.r
      # B.N.tay.expansion = B.k.min.1

      # Save for use in second derivative computation ***
      SS.obj.i$d.histories[[r]]$dB.k.squaring.history[[1]] = dB.k.min.1.r


      for(k in 2:m){ # start addition of following terms from k = 2, since k = 1 is starting point

        dB.k.r = 1/k * 1/2^N * ( SS.obj.i$B.k.squaring.history[[k-1]] %*% dA[,,r] + dB.k.min.1.r %*% A )
        # B.k = 1/k * 1/2^N * B.k.min.1 %*% A

        dB.N.tay.expansion.r = dB.N.tay.expansion.r + dB.k.r
        # B.N.tay.expansion = B.N.tay.expansion + B.k

        dB.k.min.1.r = dB.k.r
        # B.k.min.1 = B.k

        # Save for use in second derivative computation ***
        SS.obj.i$d.histories[[r]]$dB.k.squaring.history[[k]] = dB.k.min.1.r
      }

      # *** Computation of B_{,r}^[0] *** (according to recursive formula proposed in Young (2004) and Zhong & Williams (1994))
      dB.k.r = dB.N.tay.expansion.r
      # B.k = B.N.tay.expansion

      # Save for use in second derivative computation ***
      SS.obj.i$d.histories[[r]]$dB.k.scaling.history[[1]] = dB.k.r

      for(k in N:1){
        dB.k.r = 2 * dB.k.r + dB.k.r %*% SS.obj.i$B.k.scaling.history[[N-k+1]] + SS.obj.i$B.k.scaling.history[[N-k+1]] %*% dB.k.r
        # B.k = 2 * B.k + B.k %*% B.k

        # Save for use in second derivative computation ***
        SS.obj.i$d.histories[[r]]$dB.k.scaling.history[[N-k+2]] = dB.k.r
      }

      # Finally, the matrix exponential follows from formula: dB^[0] = d exp(A)
      dPmatr.i[,,r] = dB.k.r



    }




  }


  list(dPmatr = dPmatr.i,
       SS.obj = SS.obj.i)

}
