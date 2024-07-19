


dcopFun.d2 = function(m1, m2, phi, type = 'C0'){
  return( m2^(-phi-1) * (m1^(-phi) + m2^(-phi) - 1)^(-1/phi-1) )
}
