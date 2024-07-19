# # This function computes the second derivative of trans probability matrix P given
# - a transition intensity matrix Q
# - the time lag
# - the derivative of the Q matrix
# - the second derivative of the Q matrix

# For now implemented only for method = 'eigendecomp'. (Note: perhaps no need for dQ and d2Q when method is
# not eigendecom so move them as done for decomp.obj.i to the last arguments. Not sure about this yet though,
# THINK ABOUT IT).

d2P.matr.comp = function(nstates, l.params, timelag, dQmatr.i, d2Qmatr.i, start.pos.par, method = 'eigendecomp',
                         decomp.obj.i = NULL, SS.obj.i = NULL, Qmatr.i = NULL, Pmatr.i = NULL, dPmatr.i = NULL,
                         comp.par.mapping = NULL){


  if(method == 'eigendecomp'){

    tol = 12 #30 #15 # tolerance for when we are checking whether the eigenvalues are distinct... too small, not small enough? Keep an eye on this

    eigen.vec.i = decomp.obj.i$vectors
    eigen.val.i = decomp.obj.i$values
    eigen.vec.inv.i = decomp.obj.i$vectors.inv

    d2Pmatr.i = array(0, dim = c(nstates, nstates, l.params, l.params))

    # **** Check for distinct eigenvalues **** perhaps not needed since shoudl get blocked at the level of P.matr.comp() function (which is always called before)
    if( FALSE
        # !any(duplicated(round(eigen.val.i, tol)))
        ){
      # stop('The eigenvalues of the Q matrix are not distinct, the eigendecompisition method cannot be used.')

      # *** NEW! Setup of required quantities (outside of for loop) ***************
      # Computation of U.r1.r2 term (exactly the same as term in dP but with d2Q as middle matrix instead)
      # First matrix (diff of exponentials of eigenvals with diag)
      start.new = Sys.time()
      eem1 = matrix(exp(rep(eigen.val.i, nstates)*timelag), ncol = nstates)
      eem2 = t(eem1)
      diag(eem2) = 0
      diff.exp.val = eem1 - eem2

      # Need this for later
      diag(eem1) = 0
      t.exp.val.no.diag.l = timelag*eem1
      t.exp.val.no.diag.m = timelag*t(eem1)

      diff.exp.val.no.diag = eem1 - eem2
      # *****

      # Second matrix (diff of eigenvals with diag)
      eem1 = matrix(rep(eigen.val.i, nstates), ncol = nstates)
      diag(eem1) = 1/timelag
      eem2 = t(eem1)
      diag(eem2) = 0

      diff.val = eem1 - eem2

      # Need this for later
      diag(eem1) = 1 # just so don't have division times 0 later
      diff.val.no.diag = eem1 - eem2
      # ******

      # Computation of V.r1.r2 and V.r2.r1 terms (some, the rest depend on r1 and r2 and need
      # to stay within the for loop)
      frac.diff.exp.val.no.diag = diff.exp.val.no.diag/diff.val.no.diag

      frac.diff.exp.val2.no.diag = diff.exp.val.no.diag/diff.val.no.diag^2
      frac.t.exp.val.no.diag.l = t.exp.val.no.diag.l/diff.val.no.diag
      frac.t.exp.val.no.diag.m = t.exp.val.no.diag.m/diff.val.no.diag

      # diag(diff.val.no.diag) = 0 # bring it back to zero because fraction has been done # MODIFIED 05/12/2021, no need for this!! Only ever use it as denominator of division so we NEED to have some number on diagonal of this matrix (numerator always has 0 diagonal)
      # ****************************************************************************

      # To map back positions ****
      if(is.null(comp.par.mapping)) comp.par = 1 # will map back by sort of reversing the way I compacted it in the first place with comp.par if comp.par.mapping is not available
      # ***************************


      for(r1 in 1:l.params){

        G.r1 = eigen.vec.inv.i %*% dQmatr.i[,,r1] %*% eigen.vec.i

        # Define other necessary quantities
        G.r1.t = G.r1
        diag(G.r1.t) = 0

        et.m1 = matrix(rep(diag(G.r1), nstates), ncol = nstates)
        diag(et.m1) = 0  # external terms matrix 1
        # **********************************


        for(r2 in r1:l.params){ # just do upper triangle...

          G.r2 = eigen.vec.inv.i %*% dQmatr.i[,,r2] %*% eigen.vec.i

          if( is.null(comp.par.mapping) ){
            # *** COMPACT IMPLEMENTATION FIX (remember to uncomment comp.par as well) **********************
            if( max(which(r1 - start.pos.par >= 0)) == max(which(r2 - start.pos.par >= 0)) ) { # ... but remember block-diagonal structure, this checks to which transition intensity the parameter belongs, if not same for r1 and r2 then d2Q is just matrix of zeroes
              G.r1.r2 = eigen.vec.inv.i %*% d2Qmatr.i[, , comp.par] %*% eigen.vec.i
              comp.par = comp.par + 1
            } else {
              G.r1.r2 = matrix(0, nrow = nstates, ncol = nstates)
            }
            # *************************************************************************************************
          } else {

            if( !is.na(comp.par.mapping[r1, r2]) ){
              G.r1.r2 = eigen.vec.inv.i %*% d2Qmatr.i[, , comp.par.mapping[r1, r2]] %*% eigen.vec.i
            } else {
              G.r1.r2 = matrix(0, nrow = nstates, ncol = nstates)
            }

          }



          # # ******

          U.r1.r2 = G.r1.r2 * diff.exp.val/diff.val


          # Computation of V.r1.r2 and V.r2.r1 terms

          G.r2.t = G.r2
          diag(G.r2.t) = 0

          sum.terms = ( G.r1.t * frac.diff.exp.val.no.diag  ) %*% ( G.r2.t / diff.val.no.diag ) - ( G.r1.t %*% ( G.r2.t / diff.val.no.diag ) ) * frac.diff.exp.val.no.diag
          diag(sum.terms) = 0 # correction as of 02/08/2022 - this ensures that we are not contributing to the diagonal of V.r1.r2 since this is solely given by the the diag.terms term

          et.m2 = t(matrix(rep(diag(G.r2), nstates), ncol = nstates))
          diag(et.m2) = 0  # external terms matrix 2


          external.terms = et.m1 * G.r2 * (frac.t.exp.val.no.diag.l  - frac.diff.exp.val2.no.diag) + G.r1 * et.m2 * (frac.diff.exp.val2.no.diag - frac.t.exp.val.no.diag.m)

          # Now the diagonal terms # NOTE: on 11/12/2021 added additional diag() to first term, previously not present and was wrong (the first addend was a vector this way) but never noticed because always vector of 0s so far
          # diag.terms = diag( (G.r1.t * (frac.t.exp.val.no.diag.l  - frac.diff.exp.val2.no.diag) ) %*% G.r2.t ) + 0.5 * diag(diag(G.r1)) * diag(diag(G.r2)) * diag(timelag^2 * exp(eigen.val.i*timelag))
          diag.terms = diag( diag( (G.r1.t * (frac.t.exp.val.no.diag.l  - frac.diff.exp.val2.no.diag) ) %*% G.r2.t ) + 0.5 * diag(G.r1) * diag(G.r2) * timelag^2 * exp(eigen.val.i*timelag) )

          V.r1.r2 = sum.terms + external.terms + diag.terms


          # *** Now same for other way around
          sum.terms = ( G.r2.t * frac.diff.exp.val.no.diag  ) %*% ( G.r1.t / diff.val.no.diag ) - ( G.r2.t %*% ( G.r1.t / diff.val.no.diag ) ) * frac.diff.exp.val.no.diag
          diag(sum.terms) = 0 # correction as of 02/08/2022 - this ensures that we are not contributing to the diagonal of V.r2.r1 since this is solely given by the the diag.terms term

          external.terms = t(et.m2) * G.r1 * (frac.t.exp.val.no.diag.l  - frac.diff.exp.val2.no.diag) + G.r2 * t(et.m1) * (frac.diff.exp.val2.no.diag - frac.t.exp.val.no.diag.m)

          # Now the diagonal terms # NOTE: on 11/12/2021 added additional diag() to first term, previously not present and was wrong (the first addend was a vector this way) but never noticed because always vector of 0s so far
          # diag.terms = diag( (G.r2.t * (frac.t.exp.val.no.diag.l  - frac.diff.exp.val2.no.diag) ) %*% G.r1.t ) + 0.5 * diag(diag(G.r2)) * diag(diag(G.r1)) * diag(timelag^2 * exp(eigen.val.i*timelag))
          diag.terms = diag( diag( (G.r2.t * (frac.t.exp.val.no.diag.l  - frac.diff.exp.val2.no.diag) ) %*% G.r1.t ) + 0.5 * diag(G.r2) * diag(G.r1) * timelag^2 * exp(eigen.val.i*timelag) )


          V.r2.r1 = sum.terms + external.terms + diag.terms


          d2Pmatr.i[ , , r1, r2] = eigen.vec.i %*% (U.r1.r2 + V.r1.r2 + V.r2.r1) %*% eigen.vec.inv.i

        }
      }

    } else {


      E.matr = matrix(0, ncol = nstates, nrow = nstates)

      for(l in 1:nstates){ # we only need the upper triangle since its symmetric
        for(m in l:nstates){
          if(round(abs(eigen.val.i[l] - eigen.val.i[m]), tol) != 0){
            E.matr[l,m] = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/(eigen.val.i[l] - eigen.val.i[m])
          } else {
            E.matr[l,m] = timelag * exp(eigen.val.i[l]*timelag)
          }
        }
      }

      # obtain the full matrix
      E.matr = E.matr + t(E.matr)

      # add the diagonal
      diag(E.matr) = timelag * exp(eigen.val.i*timelag)

      # To map back positions ****
      if(is.null(comp.par.mapping)) comp.par = 1 # will map back by sort of reversing the way I compacted it in the first place with comp.par if comp.par.mapping is not available
      # ***************************


      for(r1 in 1:l.params){

        G.r1 = eigen.vec.inv.i %*% dQmatr.i[,,r1] %*% eigen.vec.i


        for(r2 in r1:l.params){ # just do upper triangle...

          G.r2 = eigen.vec.inv.i %*% dQmatr.i[,,r2] %*% eigen.vec.i

          if( is.null(comp.par.mapping) ){
            # *** COMPACT IMPLEMENTATION FIX (remember to uncomment comp.par as well) **********************
            if( max(which(r1 - start.pos.par >= 0)) == max(which(r2 - start.pos.par >= 0)) ) { # ... but remember block-diagonal structure, this checks to which transition intensity the parameter belongs, if not same for r1 and r2 then d2Q is just matrix of zeroes
              G.r1.r2 = eigen.vec.inv.i %*% d2Qmatr.i[, , comp.par] %*% eigen.vec.i
              comp.par = comp.par + 1
            } else {
              G.r1.r2 = matrix(0, nrow = nstates, ncol = nstates)
            }
            # *************************************************************************************************
          } else {

            if( !is.na(comp.par.mapping[r1, r2]) ){
              G.r1.r2 = eigen.vec.inv.i %*% d2Qmatr.i[, , comp.par.mapping[r1, r2]] %*% eigen.vec.i
            } else {
              G.r1.r2 = matrix(0, nrow = nstates, ncol = nstates)
            }

          }



          # # ******

          U.r1.r2 = G.r1.r2 * E.matr


          # Computation of V.r1.r2 and V.r2.r1 terms
          V.r1.r2 = matrix(0, ncol = nstates, nrow = nstates)
          V.r2.r1 = matrix(0, ncol = nstates, nrow = nstates)

          for(l in 1:nstates){
            for(m in 1:nstates){
              if(round(abs(eigen.val.i[l] - eigen.val.i[m]), tol) != 0){
                YY = 1:nstates

                for(yy in YY){
                  if(round(abs(eigen.val.i[l] - eigen.val.i[yy]), tol) != 0 & round(abs(eigen.val.i[m] - eigen.val.i[yy]), tol) != 0){
                    term1 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[yy]*timelag))/((eigen.val.i[l] - eigen.val.i[yy])*(eigen.val.i[yy] - eigen.val.i[m]))
                    term2 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/((eigen.val.i[l] - eigen.val.i[m])*(eigen.val.i[yy] - eigen.val.i[m]))

                    #print(paste(l, ',', m, ',', yy, '\n'))
                    #print(term1)
                    #print('\n')
                    #print(term2)
                    #print('\n')

                    V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * ( term1 - term2 )
                    V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * ( term1 - term2 )
                  } else if(round(abs(eigen.val.i[l] - eigen.val.i[yy]), tol) == 0 & round(abs(eigen.val.i[m] - eigen.val.i[yy]), tol) != 0){
                    term3 = timelag * exp(eigen.val.i[l]*timelag)/(eigen.val.i[l] - eigen.val.i[m])
                    term4 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/(eigen.val.i[l] - eigen.val.i[m])^2

                    #print(paste(l, ',', m, ',', yy, '\n'))
                    #print(term3)
                    #print('\n')
                    #print(term4)
                    #print('\n')

                    V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * (term3 - term4)
                    V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * (term3 - term4)
                  } else if(round(abs(eigen.val.i[l] - eigen.val.i[yy]), tol) != 0 & round(abs(eigen.val.i[m] - eigen.val.i[yy]), tol) == 0){
                    term4 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/(eigen.val.i[l] - eigen.val.i[m])^2
                    term5 = timelag * exp(eigen.val.i[m]*timelag)/(eigen.val.i[l] - eigen.val.i[m])

                    #print(paste(l, ',', m, ',', yy, '\n'))
                    #print(term4)
                    #print('\n')
                    #print(term5)
                    #print('\n')

                    V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * (term4 - term5)
                    V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * (term4 - term5)
                  }
                }


                # YY = YY[-c(l, m)]
                #
                # # V.r1.r2 term
                # for(yy in YY){
                #   term1 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[yy]*timelag))/((eigen.val.i[l] - eigen.val.i[yy])*(eigen.val.i[yy] - eigen.val.i[m]))
                #   term2 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/((eigen.val.i[l] - eigen.val.i[m])*(eigen.val.i[yy] - eigen.val.i[m]))
                #   V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * ( term1 - term2 )
                # }
                # term3 = timelag * exp(eigen.val.i[l]*timelag)/(eigen.val.i[l] - eigen.val.i[m])
                # term4 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/(eigen.val.i[l] - eigen.val.i[m])^2
                # term5 = timelag * exp(eigen.val.i[m]*timelag)/(eigen.val.i[l] - eigen.val.i[m])
                # V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,l] * G.r2[l,m] * ( term3 - term4 ) + G.r1[l,m] * G.r2[m,m] * (term4 - term5)
                #
                # # V.r2.r1 term
                # for(yy in YY){
                #   term1 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[yy]*timelag))/((eigen.val.i[l] - eigen.val.i[yy])*(eigen.val.i[yy] - eigen.val.i[m]))
                #   term2 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[m]*timelag))/((eigen.val.i[l] - eigen.val.i[m])*(eigen.val.i[yy] - eigen.val.i[m]))
                #   V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * (term1 - term2)
                # }
                # V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,l] * G.r1[l,m] * (term3 - term4) + G.r2[l,m] * G.r1[m,m] * (term4 - term5)
              } else {
                YY = 1:nstates

                for(yy in YY){
                  if(round(abs(eigen.val.i[l] - eigen.val.i[yy]), tol) != 0 & round(abs(eigen.val.i[m] - eigen.val.i[yy]), tol) != 0){
                    term6 = timelag * exp(eigen.val.i[l]*timelag)/(eigen.val.i[l] - eigen.val.i[yy])
                    term7 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[yy]*timelag))/(eigen.val.i[l] - eigen.val.i[yy])^2

                    #print(paste(l, ',', m, ',', yy, '\n'))
                    #print(term6)
                    #print('\n')
                    #print(term7)
                    #print('\n')

                    V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * ( term6 - term7 )
                    V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * ( term6 - term7 )
                  } else if(round(abs(eigen.val.i[l] - eigen.val.i[yy]), tol) == 0 & round(abs(eigen.val.i[m] - eigen.val.i[yy]), tol) == 0){
                    term8 = 0.5 * timelag^2 * exp(eigen.val.i[yy]*timelag)

                    #print(paste(l, ',', m, ',', yy, '\n'))
                    #print(term8)
                    #print('\n')


                    V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * term8
                    V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * term8
                  }
                }

                # YY = YY[-unique(c(l, m))]
                #
                # # V.r1.r2 term
                # for(yy in YY){
                #   term6 = timelag * exp(eigen.val.i[l]*timelag)/(eigen.val.i[l] - eigen.val.i[yy])
                #   term7 = (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[yy]*timelag))/(eigen.val.i[l] - eigen.val.i[yy])^2
                # } V.r1.r2[l,m] = V.r1.r2[l,m] + G.r1[l,yy] * G.r2[yy,m] * ( term6 - term7 )
                # for(yy in unique(c(l,m))) V.r1.r2[l,m] = V.r1.r2[l,m] + 0.5 * G.r1[l,yy] * G.r2[yy,m] * timelag^2 * exp(eigen.val.i[yy]*timelag)
                #
                # # V.r2.r1 term
                # for(yy in YY) V.r2.r1[l,m] = V.r2.r1[l,m] + G.r2[l,yy] * G.r1[yy,m] * ( timelag * exp(eigen.val.i[l]*timelag)/(eigen.val.i[l] - eigen.val.i[yy]) - (exp(eigen.val.i[l]*timelag) - exp(eigen.val.i[yy]*timelag))/(eigen.val.i[l] - eigen.val.i[yy])^2 )
                # for(yy in unique(c(l,m))) V.r2.r1[l,m] = V.r2.r1[l,m] + 0.5 * G.r2[l,yy] * G.r1[yy,m] * timelag^2 * exp(eigen.val.i[yy]*timelag)
              }
            }
          }


          d2Pmatr.i[ , , r1, r2] = eigen.vec.i %*% (U.r1.r2 + V.r1.r2 + V.r2.r1) %*% eigen.vec.inv.i

        }
      }


    }




  }


  if( method == 'analytic'){

    # allow user to define this tolerance? Does it create overflow issues? When?
    tol = 1e-15 # how small a difference counts as equality?

    d2Pmatr.i = array(0, dim = c(nstates, nstates, l.params, l.params))

    # # *** OLD NAIVE IMPLEMENTATION *******************************************************
    # # Repeated quantities ***
    # rq1 = exp(-Qmatr.i[2,3]*timelag) - exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)
    # drq1 = - timelag * dQmatr.i[2,3,] * exp(-Qmatr.i[2,3]*timelag) + timelag * (dQmatr.i[1,2,] + dQmatr.i[1,3,]) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)
    # d2rq1 = - timelag * d2Qmatr.i[2,3,,] * exp(-Qmatr.i[2,3]*timelag) + timelag^2 * ( dQmatr.i[2,3,] %*% t(dQmatr.i[2,3,]) ) * exp(-Qmatr.i[2,3]*timelag) +
    #   timelag * (d2Qmatr.i[1,2,,] + d2Qmatr.i[1,3,,]) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag) - timelag^2 * ( (dQmatr.i[1,2,] + dQmatr.i[1,3,]) %*% t(dQmatr.i[1,2,] + dQmatr.i[1,3,]) ) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)
    #
    # rq2 = Qmatr.i[1,2] + Qmatr.i[1,3] - Qmatr.i[2,3]
    # drq2 = dQmatr.i[1,2,] + dQmatr.i[1,3,] - dQmatr.i[2,3,]
    # d2rq2 = d2Qmatr.i[1,2,,] + d2Qmatr.i[1,3,,] - d2Qmatr.i[2,3,,]
    # # *****
    #
    #
    # d2Pmatr.i[1,1,,] = timelag * ( dPmatr.i[1,1,]  %*% t(dQmatr.i[1,1,]) ) + Pmatr.i[1,1] * timelag * d2Qmatr.i[1,1,,]
    #
    # if( abs( (Qmatr.i[1,2] + Qmatr.i[1,3]) - Qmatr.i[2,3] ) > tol ){
    #   d2Pmatr.i[1,2,,] = ( rq2 * ( d2Qmatr.i[1,2,,] * rq1 + ( dQmatr.i[1,2,] %*% t(drq1) ) ) - rq1 * ( dQmatr.i[1,2,] %*% t(drq2) ) ) / rq2^2 -
    #     ( rq2^2 * ( rq1 * ( ( dQmatr.i[1,2,] %*% t(drq2) ) + Qmatr.i[1,2] * d2rq2 ) + Qmatr.i[1,2] * ( drq1 %*% t(drq2) ) ) - 2 * Qmatr.i[1,2] * rq1 * rq2 * ( drq2 %*% t(drq2) ) ) / rq2^4 +
    #     ( rq2 * ( dQmatr.i[1,2,] %*% t(drq1) ) - Qmatr.i[1,2] * ( drq2 %*% t(drq1) ) ) / rq2^2 + Qmatr.i[1,2]/rq2 * d2rq1
    # } else {
    #   d2Pmatr.i[1,2,,] = - timelag * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag) * ( (dQmatr.i[1,2,] + dQmatr.i[1,3,]) %*% t( timelag * dQmatr.i[1,2,] - timelag^2 * Qmatr.i[1,2] * (dQmatr.i[1,2,] + dQmatr.i[1,3,]) ) ) +
    #     exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag) * ( timelag * d2Qmatr.i[1,2,,]  - timelag^2 * ( dQmatr.i[1,2,] %*% t(dQmatr.i[1,2,] + dQmatr.i[1,3,]) ) - timelag^2 * Qmatr.i[1,2] * (d2Qmatr.i[1,2,,] + d2Qmatr.i[1,3,,]) )
    # }
    #
    # d2Pmatr.i[1,3,,] = - (d2Pmatr.i[1,1,,] + d2Pmatr.i[1,2,,])
    #
    # d2Pmatr.i[2,1,,] = 0
    # d2Pmatr.i[2,2,,] = - timelag * dPmatr.i[2,2,] %*% t(dQmatr.i[2,3,]) - timelag * Pmatr.i[2,2] * d2Qmatr.i[2,3,,]
    # d2Pmatr.i[2,3,,] = - (d2Pmatr.i[2,1,,] + d2Pmatr.i[2,2,,])
    #
    # d2Pmatr.i[3,1,,] = 0
    # d2Pmatr.i[3,2,,] = 0
    # d2Pmatr.i[3,3,,] = 0
    # # *************************************************************************************************************************


    # **** NEW IMPLEMENTATION FOLLOWING FROM UPDATED FORM GIVEN TO d2Q (no longer nstates x nstates x l.params x l.params x nrow(data) but compact version)
    # In practice: we will bring it back to that form.

    d2Qmatr.i.comp = d2Qmatr.i # save compact one differently
    d2Qmatr.i = array(0, dim = c(nstates, nstates, l.params, l.params))

    # To map back positions *******************************************************************************
    if(is.null(comp.par.mapping)) comp.par = 1 # will map back by sort of reversing the way I compacted it in the first place with comp.par if comp.par.mapping is not available
    # ***************************

    if( is.null(comp.par.mapping) ){
      for(r1 in 1:l.params){
        for(r2 in r1:l.params){
      # *** COMPACT IMPLEMENTATION FIX (remember to uncomment comp.par as well) **********************
      if( max(which(r1 - start.pos.par >= 0)) == max(which(r2 - start.pos.par >= 0)) ) { # ... but remember block-diagonal structure, this checks to which transition intensity the parameter belongs, if not same for r1 and r2 then d2Q is just matrix of zeroes
        d2Qmatr.i[, , r1, r2] = d2Qmatr.i.comp[, , comp.par]
        comp.par = comp.par + 1
      } else {
        d2Qmatr.i[, , r1, r2] = matrix(0, nrow = nstates, ncol = nstates)
      }
      # *************************************************************************************************
        }}
    } else {
      for(r1 in 1:l.params){
        for(r2 in r1:l.params){
      if( !is.na(comp.par.mapping[r1, r2]) ){
        d2Qmatr.i[, , r1, r2] = d2Qmatr.i.comp[, , comp.par.mapping[r1, r2]]
      } else {
        d2Qmatr.i[, , r1, r2] = matrix(0, nrow = nstates, ncol = nstates)
      }
        }}
    }




    # ******************************************************************************************************


    # Repeated quantities ***
    rq1 = exp(-Qmatr.i[2,3]*timelag) - exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)
    drq1 = - timelag * dQmatr.i[2,3,] * exp(-Qmatr.i[2,3]*timelag) + timelag * (dQmatr.i[1,2,] + dQmatr.i[1,3,]) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)
    d2rq1 = - timelag * d2Qmatr.i[2,3,,] * exp(-Qmatr.i[2,3]*timelag) + timelag^2 * ( dQmatr.i[2,3,] %*% t(dQmatr.i[2,3,]) ) * exp(-Qmatr.i[2,3]*timelag) +
      timelag * (d2Qmatr.i[1,2,,] + d2Qmatr.i[1,3,,]) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag) - timelag^2 * ( (dQmatr.i[1,2,] + dQmatr.i[1,3,]) %*% t(dQmatr.i[1,2,] + dQmatr.i[1,3,]) ) * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag)

    rq2 = Qmatr.i[1,2] + Qmatr.i[1,3] - Qmatr.i[2,3]
    drq2 = dQmatr.i[1,2,] + dQmatr.i[1,3,] - dQmatr.i[2,3,]
    d2rq2 = d2Qmatr.i[1,2,,] + d2Qmatr.i[1,3,,] - d2Qmatr.i[2,3,,]
    # *****


    d2Pmatr.i[1,1,,] = timelag * ( dPmatr.i[1,1,]  %*% t(dQmatr.i[1,1,]) ) + Pmatr.i[1,1] * timelag * d2Qmatr.i[1,1,,]
    d2Pmatr.i[1,1,,][lower.tri(d2Pmatr.i[1,1,,], diag = FALSE)] = 0

    if( abs( (Qmatr.i[1,2] + Qmatr.i[1,3]) - Qmatr.i[2,3] ) > tol ){
      d2Pmatr.i[1,2,,] = ( rq2 * ( d2Qmatr.i[1,2,,] * rq1 + ( drq1 %*% t(dQmatr.i[1,2,]) ) ) - rq1 * ( drq2 %*% t(dQmatr.i[1,2,]) ) ) / rq2^2 -
        ( rq2^2 * ( rq1 * ( ( dQmatr.i[1,2,] %*% t(drq2) ) + Qmatr.i[1,2] * d2rq2 ) + Qmatr.i[1,2] * ( drq1 %*% t(drq2) ) ) - 2 * Qmatr.i[1,2] * rq1 * rq2 * ( drq2 %*% t(drq2) ) ) / rq2^4 +
        ( rq2 * ( dQmatr.i[1,2,] %*% t(drq1) ) - Qmatr.i[1,2] * ( drq2 %*% t(drq1) ) ) / rq2^2 + Qmatr.i[1,2]/rq2 * d2rq1
    } else {
      d2Pmatr.i[1,2,,] = - timelag * exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag) * ( (dQmatr.i[1,2,] + dQmatr.i[1,3,]) %*% t( timelag * dQmatr.i[1,2,] - timelag^2 * Qmatr.i[1,2] * (dQmatr.i[1,2,] + dQmatr.i[1,3,]) ) ) +
        exp(-(Qmatr.i[1,2] + Qmatr.i[1,3])*timelag) * ( timelag * d2Qmatr.i[1,2,,]  - timelag^2 * ( dQmatr.i[1,2,] %*% t(dQmatr.i[1,2,] + dQmatr.i[1,3,]) ) - timelag^2 * Qmatr.i[1,2] * (d2Qmatr.i[1,2,,] + d2Qmatr.i[1,3,,]) )
    }
    d2Pmatr.i[1,2,,][lower.tri(d2Pmatr.i[1,2,,], diag = FALSE)] = 0


    d2Pmatr.i[1,3,,] = - (d2Pmatr.i[1,1,,] + d2Pmatr.i[1,2,,]) # already cleaned lower triangle

    d2Pmatr.i[2,1,,] = 0
    d2Pmatr.i[2,2,,] = - timelag * dPmatr.i[2,2,] %*% t(dQmatr.i[2,3,]) - timelag * Pmatr.i[2,2] * d2Qmatr.i[2,3,,]
    d2Pmatr.i[2,2,,][lower.tri(d2Pmatr.i[2,2,,], diag = FALSE)] = 0

    d2Pmatr.i[2,3,,] = - (d2Pmatr.i[2,1,,] + d2Pmatr.i[2,2,,])

    d2Pmatr.i[3,1,,] = 0
    d2Pmatr.i[3,2,,] = 0
    d2Pmatr.i[3,3,,] = 0

    # ***************************************************************************************************************************************************
  }



  if( method == 'scaling&squaring'){

    warning('The scaling&squaring method is still a work in progress. Please use the other methods instead as they are faster.')


    # Values of N and m proposed by Zhong & Williams (1994)
    N = 20
    m = 4

    # Matrix setup
    A = Qmatr.i*timelag
    dA = dQmatr.i*timelag
    d2A = d2Qmatr.i*timelag
    Id.m = diag(dim(A)[1])

    if( 1/2^N * max(abs(A)) >= 1  ){
      N = ceiling(log2(max(abs(A))))
      write(paste('Largest Q value found was', max(abs(A)), 'so N was set to', N, '\n'), file = 'testing_adaptive_S&S.txt', append = TRUE)
    }


    d2Pmatr.i = array(0, dim = c(nstates, nstates, l.params, l.params))

    # To map back positions ****
    comp.par = 1 # will map back by sort of reversing the way I compacted it in the first place
    # ***************************

    for(r1 in 1:l.params){
      for(r2 in r1:l.params){ # just do upper triangle...


        # *** COMPACT IMPLEMENTATION FIX (remember to uncomment comp.par as well) **********************
        if( max(which(r1 - start.pos.par >= 0)) == max(which(r2 - start.pos.par >= 0)) ) { # ... but remember block-diagonal structure, this checks to which transition intensity the parameter belongs, if not same for r1 and r2 then d2Q is just matrix of zeroes
          d2A.r1.r2 = d2A[, , comp.par]
          comp.par = comp.par + 1
        } else {
          d2A.r1.r2 = matrix(0, nrow = nstates, ncol = nstates)
        }
        # *************************************************************************************************


        # *** Computation of B_{,r1 r2}^[N] *** as proposed in Young (2004), with correction of typos
        # B.k.min.1 = 1/2^N * A
        # dB.k.min.1.r1 = 1/2^N * dA[,,r1] # k = 1 term of summation
        # dB.k.min.1.r2 = 1/2^N * dA[,,r2] # k = 1 term of summation
        d2B.k.min.1.r1.r2 = 1/2^N * d2A.r1.r2

        d2B.N.tay.expansion.r1.r2 = d2B.k.min.1.r1.r2
        # dB.N.tay.expansion.r1 = dB.k.min.1.r1
        # dB.N.tay.expansion.r2 = dB.k.min.1.r2
        # B.N.tay.expansion = B.k.min.1

        for(k in 2:m){ # start addition of following terms from k = 2, since k = 1 is starting point

          d2B.k.r1.r2 = 1/k * 1/2^N *( SS.obj.i$B.k.squaring.history[[k-1]] %*% d2A.r1.r2 + SS.obj.i$d.histories[[r1]]$dB.k.squaring.history[[k-1]] %*% dA[,,r2] +
                                         SS.obj.i$d.histories[[r2]]$dB.k.squaring.history[[k-1]] %*% dA[,,r1] + d2B.k.min.1.r1.r2 %*% A )
          # dB.k.r1 = 1/k * 1/2^N * ( B.k.min.1 %*% dA[,,r1] + dB.k.min.1.r1 %*% A )
          # dB.k.r2 = 1/k * 1/2^N * ( B.k.min.1 %*% dA[,,r2] + dB.k.min.1.r2 %*% A )
          # B.k = 1/k * 1/2^N * B.k.min.1 %*% A

          d2B.N.tay.expansion.r1.r2 = d2B.N.tay.expansion.r1.r2 + d2B.k.r1.r2
          # dB.N.tay.expansion.r1 = dB.N.tay.expansion.r1 + dB.k.r1
          # dB.N.tay.expansion.r2 = dB.N.tay.expansion.r2 + dB.k.r2
          # B.N.tay.expansion = B.N.tay.expansion + B.k

          d2B.k.min.1.r1.r2 = d2B.k.r1.r2
          # dB.k.min.1.r1 = dB.k.r1
          # dB.k.min.1.r2 = dB.k.r2
          # B.k.min.1 = B.k
        }

        # *** Computation of B_{,r1 r2}^[0] *** (according to recursive formula proposed in Young (2004) and Zhong & Williams (1994))
        d2B.k.r1.r2 = d2B.N.tay.expansion.r1.r2
        # dB.k.r1 = dB.N.tay.expansion.r1
        # dB.k.r2 = dB.N.tay.expansion.r2
        # B.k = B.N.tay.expansion

        for(k in N:1){
          # Use all recycled quantities
          d2B.k.r1.r2 = 2 * d2B.k.r1.r2 + d2B.k.r1.r2 %*% SS.obj.i$B.k.scaling.history[[N-k+1]] +
            SS.obj.i$d.histories[[r1]]$dB.k.scaling.history[[N-k+1]] %*% SS.obj.i$d.histories[[r2]]$dB.k.scaling.history[[N-k+1]] +
            SS.obj.i$d.histories[[r2]]$dB.k.scaling.history[[N-k+1]] %*% SS.obj.i$d.histories[[r1]]$dB.k.scaling.history[[N-k+1]] +
            SS.obj.i$B.k.scaling.history[[N-k+1]] %*% d2B.k.r1.r2
          # dB.k.r1 = 2 * dB.k.r1 + dB.k.r1 %*% B.k + B.k %*% dB.k.r1
          # dB.k.r2 = 2 * dB.k.r2 + dB.k.r2 %*% B.k + B.k %*% dB.k.r2
          # B.k = 2 * B.k + B.k %*% B.k
        }

        # Finally, the matrix exponential follows from formula: B^[0] + I = exp(A)
        d2Pmatr.i[,,r1,r2] = d2B.k.r1.r2


      }
    }


  }

  d2Pmatr.i


}



