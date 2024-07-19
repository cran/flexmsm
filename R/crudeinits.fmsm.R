

# adaptation of the function msm::crudeinits.msm()

crudeinits.fmsm = function(state, time, subject, qmatrix, data = NULL)
{

  # mf <- model.frame(formula, data = data, na.action = NULL)
  state <- data[,state] #mf[,1]
  # if (is.factor(state))
  #   state <- as.numeric(as.character(state))
  time <- data[,time] # mf[, 2]
  n <- length(state)
  if (missing(subject))
    subject <- rep(1, n)
  else if (!is.null(data))
    subject <- data[,subject]
  if (is.null(subject))
    subject <- rep(1, n)
  notna <- !is.na(subject) & !is.na(time) & !is.na(state)
  subject <- subject[notna]
  time <- time[notna]
  state <- state[notna]


  n <- length(state)
  lastsubj <- !duplicated(subject, fromLast = TRUE)
  timecontrib <- ifelse(lastsubj, NA, c(time[2:n], 0) - time)
  tottime <- tapply(timecontrib[!lastsubj], state[!lastsubj],
                    sum)

  ntrans <- statetable.fmsm(state, subject, data = NULL)
  nst <- nrow(qmatrix)
  estmat <- matrix(0, nst, nst)
  rownames(estmat) <- colnames(estmat) <- paste(1:nst)
  tab <- sweep(ntrans, 1, tottime, "/")
  for (i in 1:nst) for (j in 1:nst) if ((paste(i) %in% rownames(tab)) &&
                                        (paste(j) %in% colnames(tab)))
    estmat[paste(i), paste(j)] <- tab[paste(i), paste(j)]
  estmat[qmatrix == 0] <- 0

  diag(estmat) <- 0
  diag(estmat) <- -rowSums(estmat)

  rownames(estmat) <- rownames(qmatrix)
  colnames(estmat) <- colnames(qmatrix)

  estmat
}
