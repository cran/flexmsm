#' Flexible transition intensity based models for univariate multistate processes
#'
#' @param x Fitted model object.
#' @param digits Number of digits printed in the output.
#' @param signif.stars By default significance stars are printed alongside the output.
#' @param ... Other arguments.
#'
#'
#' @return Prints model term summaries.
#' @export
#'
#'
print.summary.fmsm = function(x, digits = max(3, getOption("digits") - 3),
                             signif.stars = getOption("show.signif.stars"), ...){


  suStf = x$suStf
  state.pairs.CT.out = x$state.pairs.CT.out

  cat('\nFlexible multistate survival model with', sum(suStf$whereQ != 0), 'transitions.\n\n')


  for(jj in 1:sum(suStf$whereQ != 0)){

    if(sum(suStf$whereQ != 0) > 15) warning('The current code only displays results for up to 15 transitions. Contact maintainer for further details.')


    cat('\nTRANSITION', x$idx.lab[jj,2], '->', x$idx.lab[jj,1], '\n')
    cat('\nFormula: ')
    print(x$short.formula[[jj]])
    cat('\n')

    cat("Parametric coefficients:\n")
    printCoefmat(x$table[[jj]], digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat('\n')

    if(!is.null(x$tableN[[jj]])) {
      cat("Smooth components' approximate significance:\n")
      printCoefmat(x$tableN[[jj]], digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
      cat('\n')

    }


    cat('n(', x$idx.lab[jj,2], '->', x$idx.lab[jj,1], ') = ', state.pairs.CT.out$counts[x$idx.lab[jj,2], x$idx.lab[jj,1]],
        '   total edf = ', round(sum(x$postpost[[jj]]$edf), 3),
        '\n\n\n', sep = '')

  }


  cat('\nn =', x$n, '  N =', x$N, '  Total edf =', round(x$t.edf,3), '\n\n')






}
