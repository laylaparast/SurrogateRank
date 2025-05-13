est.power = function(n.total, rho=0.80, u.y.alt,delta.alt, power.want.s = 0.7){
  delta.se.est = sqrt(2*(1-rho)/(3*n.total))
  sd.null = sqrt((n.total+1)/(3*n.total^2))
  z.alpha.2 = qnorm(0.975)
  u.s.power = 1/2-(qnorm(1-power.want.s)-z.alpha.2)*(sd.null)
  z.alpha = qnorm(0.95)
  est.power = pnorm((1/delta.se.est)*(u.y.alt - u.s.power-delta.alt) - z.alpha)
  return(est.power)
}