


d2copFun.dphidphi = function(m1, m2, phi, type = 'c0'){ # CHECKED FOR CORRECTNESS
  trm = m1^-phi + m2^-phi - 1
  return( trm^(-1/phi) * ( phi^-2 * log(trm) + (phi * trm)^-1 * (log(m1) * m1^-phi + log(m2) * m2^-phi) )^2
          - trm^(-1/phi) * (2 * phi^-3 * log(trm) + 2 * phi^-2 * trm^-1 * (log(m1) * m1^-phi + log(m2) * m2^-phi)
                            + (phi * trm)^-1 * (log(m1)^2 * m1^-phi + log(m2)^2 * m2^-phi) - phi^-1 * trm^-2 * (log(m1) * m1^-phi + log(m2) * m2^-phi)^2 )
  )
}
