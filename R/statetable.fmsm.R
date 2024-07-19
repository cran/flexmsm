

statetable.fmsm = function(state, subject, data = NULL){

  if (!is.null(data)) {
    data  <- as.data.frame(data)
    state <- data[,state]
  }

  n <- length(state)

  if (!is.null(data)){
    if (missing(subject)) {
      subject <-  rep(1, n)
    } else {
      subject <- data[,subject]
    }
  }

  subject <- match(subject, unique(subject))
  prevsubj <- c(NA, subject[1:(n - 1)])
  previous <- c(NA, state[1:(n - 1)])
  previous[prevsubj != subject] <- NA
  ntrans <- table(previous, state)
  names(dimnames(ntrans)) <- c("from", "to")

  ntrans

}
