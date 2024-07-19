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
print.summary.fMmsm = function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...){


  cat("\nCopula based model for two flexible dependant multistate survival processes.\n\n")


  # PROCESS 1
  suStf1 = x$suStf1
  state.pairs.CT.out1 = x$state.pairs.CT.out1

  cat("\nPROCESS 1\n")
  cat("\nFlexible multistate survival model with", sum(suStf1$whereQ != 0), "transitions.\n\n")


  ii = 1
  for (i in 1:sum(suStf1$whereQ != 0)) {
    if (sum(suStf1$whereQ != 0) > 15) warning("The current code only displays results for up to 15 transitions. Contact maintainer for further details.")
    cat("\nTRANSITION", x$idx.lab1[i, 2], "->", x$idx.lab1[i,1], "\n")
    cat("\nFormula: ")
    print(x$short.formula1[[i]])
    cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$table[[ii]], digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat("\n")

    if (!is.null(x$tableN[[ii]])) {
      cat("Smooth components' approximate significance:\n")
      printCoefmat(x$tableN[[ii]], digits = digits, signif.stars = signif.stars,na.print = "NA", ...)
      cat("\n")
    }

    cat("n(", x$idx.lab1[i, 2], "->", x$idx.lab1[i, 1], ") = ",
        state.pairs.CT.out1$counts[x$idx.lab1[i,2], x$idx.lab1[i, 1]], "   total edf = ", round(sum(x$postpost[[ii]]$edf),3), "\n\n\n", sep = "")

    ii = ii + 1
  }

  cat("\nn =", x$n1, "  Total edf =", round(sum(x$t.edf1), 3), "\n\n")


  # PROCESS 2
  suStf2 = x$suStf2
  state.pairs.CT.out2 = x$state.pairs.CT.out2

  cat("\nPROCESS 2\n")
  cat("\nFlexible multistate survival model with", sum(suStf2$whereQ != 0), "transitions.\n\n")

  for (i in 1:sum(suStf2$whereQ != 0)) {
    if (sum(suStf2$whereQ != 0) > 15) warning("The current code only displays results for up to 15 transitions. Contact maintainer for further details.")
    cat("\nTRANSITION", x$idx.lab2[i, 2], "->", x$idx.lab2[i,1], "\n")
    cat("\nFormula: ")
    print(x$short.formula2[[i]])
    cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$table[[ii]], digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat("\n")

    if (!is.null(x$tableN[[ii]])) {
      cat("Smooth components' approximate significance:\n")
      printCoefmat(x$tableN[[ii]], digits = digits, signif.stars = signif.stars,na.print = "NA", ...)
      cat("\n")
    }

    cat("n(", x$idx.lab2[i, 2], "->", x$idx.lab2[i, 1], ") = ",
        state.pairs.CT.out2$counts[x$idx.lab2[i,2], x$idx.lab2[i, 1]], "   total edf = ", round(sum(x$postpost[[ii]]$edf),3), "\n\n\n", sep = "")
  }

  cat("\nn =", x$n2, "  Total edf =", round(sum(x$t.edf2), 3), "\n\n")




  cat("\nN =", x$N, "  phi =", round(x$phiHat,4), "Total edf =", round(x$t.edf, 3), "\n\n") # add also total total edf ?? check if this is the sum of the two total edfs + 1 (for the phi) ... ?


}
