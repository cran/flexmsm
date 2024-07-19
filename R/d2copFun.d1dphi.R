


d2copFun.d1dphi = function(m1, m2, phi, type = 'c0'){ # CHECKED FOR CORRECTNESS
  return( m1^(-phi-1) * ( m1^-phi + m2^-phi - 1 )^(-1/phi-1) * ( phi^-2 * log(m1^-phi + m2^-phi - 1) - log(m1) + (phi^-1 + 1) * (log(m1) * m1^-phi + log(m2) * m2^-phi)/(m1^-phi + m2^-phi - 1) ) )
}

