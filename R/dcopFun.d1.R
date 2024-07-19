



dcopFun.d1 = function(m1, m2, phi, type = 'C0'){
  return( m1^(-phi-1) * (m1^(-phi) + m2^(-phi) - 1)^(-1/phi-1) )
}
