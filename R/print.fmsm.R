
#' Print a fmsm object
#'
#' @description
#' The print method for the \code{fmsmObject} produced by \code{\link{fmsm}}.
#'
#'
#' @param x \code{fmsm} object produced by function \code{\link{fmsm}}.
#' @param ... Unused in this case.
#'
#' @return \code{print.fmsm} prints out a matrix summarising the positions of the transition intensities, the transition intensities formulae, the total number of observations, etc for the fitted multistate survival model.
#'
#' @export print.fmsm
#' @export
#'


print.fmsm = function(x, ...){


  suStf = x$suStf

  idx.lab = which(t(suStf$whereQ) != 0, arr.ind = T)

  state.pairs.CT.out = state.pairs.CT(data = suStf$data, whereQ = suStf$whereQ, nstates = suStf$nstates,
                                      time = suStf$tte, state = '(state)', id = '(id)')

  cat('\nFlexible multistate survival model with', sum(suStf$whereQ != 0), 'transitions.\n\n')

  # cat('\nPositions of the transition intensities in the Q matrix:\n')
  # print(suStf$whereQ)
  # cat('\n')

  for(jj in 1:sum(suStf$whereQ != 0)){

    if(sum(suStf$whereQ != 0) > 15) warning('The current code only displays results for up to 15 transitions. Contact maintainer for further details.')


    cat('\nTRANSITION ', idx.lab[jj,2], ' -> ', idx.lab[jj,1], '\n', sep = '')
    cat('\nFormula: ')
    print(x$short.formula[[jj]])
    cat('\n')

    cat('n(', idx.lab[jj,2], '->', idx.lab[jj,1], ') = ', state.pairs.CT.out$counts[idx.lab[jj,2], idx.lab[jj,1]],
        '   total edf = ', round(sum(x$msm.post.object$msm.post.object[[jj]]$edf), 3),
        '\n\n\n', sep = '')



  }


  cat('\nn =', x$n, '  N =', x$N, '  Total edf =', round(x$t.edf,3), '\n\n')


}
