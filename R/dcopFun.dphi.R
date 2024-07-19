


dcopFun.dphi = function(m1, m2, phi, type = 'C0'){
  # return( ( m1^(-phi) + m2^(-phi) - 1 )^(-1/phi)
  # * ( phi^(-2) * log(m1^(-phi) + m2^(-phi) - 1) + phi^(-1) * (m1^(-phi) + m2^(-phi) - 1)^(-1) * (log(m1) * m1^(-phi) + log(m2) * m2^(-phi)) ) )
  return(((-1 + m1^-phi + m2^-phi)^(-1/phi) *((phi *(m2^phi *log(m1) + m1^phi * log(m2)))/(
    m2^phi - m1^phi* (-1 + m2^phi)) +
      log(-1 + m1^-phi + m2^-phi)) )/phi^2)
}
