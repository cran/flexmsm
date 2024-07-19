



d2copFun.d1d2 = function(m1, m2, phi, type = 'c0'){ # CHECKED FOR CORRECTNESS
  return( (1 + phi) * m1^(-phi-1) * m2^(-phi-1) * ( m1^-phi + m2^-phi - 1 )^(-1/phi-2) )
}
