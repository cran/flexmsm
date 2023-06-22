
#' Function to plot the smooths included in the model specifications.
#'
#' @param x Fitted model object.
#' @param ... Other graphical arguments.
#'
#' @return Plots the smooths.
#'
#' @export plot.fmsm
#' @export
#'


plot.fmsm = function(x, ...){


  # sub.edf <- function(lab, edf) {
  #   pos <- regexpr(":", lab)[1]
  #   if (pos < 0) {
  #     pos <- nchar(lab) - 1
  #     lab <- paste(substr(lab, start = 1, stop = pos),
  #                  ",", round(edf, digits = 2), ")", sep = "")
  #   }
  #   else {
  #     lab1 <- substr(lab, start = 1, stop = pos - 2)
  #     lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
  #     lab <- paste(lab1, ",", round(edf, digits = 2),
  #                  lab2, sep = "")
  #   }
  #   lab
  # }

  if(is.na(match('rug', names(match.call())))) rug = TRUE

  suStf = x$suStf
  idx.lab = which(t(suStf$whereQ) !=0, arr.ind = T)

  for(ii in 1:length(x$msm.post.object$msm.post.object)){

    for(jj in 1:length(x$msm.post.object$msm.post.object[[ii]]$smooth)){

      mgcv::plot.gam(x$msm.post.object$msm.post.object[[ii]],
               select = jj, rug = FALSE,
               ...)

      if(rug & length(x$msm.post.object$msm.post.object[[ii]]$smooth[[jj]]$term) == 1){ # i.e. rugplots only for 1d smooths
        pair.states = suStf$data[-nrow(suStf$data), '(state)'] == idx.lab[ii, 2] & suStf$data[-1, '(state)'] == idx.lab[ii, 1]
        same.person = suStf$data[-nrow(suStf$data), '(id)'] == suStf$data[-1, '(id)']

        act.cov = suStf$data[-1, x$msm.post.object$msm.post.object[[ii]]$smooth[[jj]]$term][pair.states & same.person]

        rug(act.cov)
      }

    }






  }




}
