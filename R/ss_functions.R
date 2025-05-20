delta.calculate <- function(full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, szero = NULL) {
  if (!is.null(full.data)) {
    yone <- full.data[full.data[, 3] == 1, 1]
    yzero <- full.data[full.data[, 3] == 0, 1]
    sone <- full.data[full.data[, 3] == 1, 2]
    szero <- full.data[full.data[, 3] == 0, 2]
  }
  test.y <- wilcox.test(yone, yzero, exact = F)
  test.s <- wilcox.test(sone, szero, exact = F)
  n1.f <- length(yone)
  n0.f <- length(yzero)
  u.y <- (n1.f * n0.f)^(-1) * test.y$statistic
  u.s <- (n1.f * n0.f)^(-1) * test.s$statistic
  delta.estimate <- u.y - u.s

  # variance, X is treated, Y is control, n is size of control, m is size of treatment
  # need to calculate S_01 and S_10
  # S_10
  # 1,1 ELEMENT
  m.count <- length(yone)
  n.count <- length(yzero)
  # need V10.Xi.Y
  V10.Xi.Y <- sapply(yone, var.wil, b = yzero)
  V10.Xi.S <- sapply(sone, var.wil, b = szero)
  V01.Yj.Y <- sapply(yzero, var.wil, b = yone, flip = TRUE)
  V01.Yj.S <- sapply(szero, var.wil, b = sone, flip = TRUE)


  s10.11.YY <- 1 / (m.count - 1) * sum((V10.Xi.Y - u.y) * (V10.Xi.Y - u.y))
  s10.12.YS <- 1 / (m.count - 1) * sum((V10.Xi.Y - u.y) * (V10.Xi.S - u.s))
  s10.22.SS <- 1 / (m.count - 1) * sum((V10.Xi.S - u.s) * (V10.Xi.S - u.s))
  s10.21.SY <- 1 / (m.count - 1) * sum((V10.Xi.Y - u.y) * (V10.Xi.S - u.s))

  s01.11.YY <- 1 / (n.count - 1) * sum((V01.Yj.Y - u.y) * (V01.Yj.Y - u.y))
  s01.12.YS <- 1 / (n.count - 1) * sum((V01.Yj.Y - u.y) * (V01.Yj.S - u.s))
  s01.22.SS <- 1 / (n.count - 1) * sum((V01.Yj.S - u.s) * (V01.Yj.S - u.s))
  s01.21.SY <- 1 / (n.count - 1) * sum((V01.Yj.Y - u.y) * (V01.Yj.S - u.s))

  S10 <- matrix(c(s10.11.YY, s10.12.YS, s10.21.SY, s10.22.SS), nrow = 2, ncol = 2, byrow = TRUE)
  S01 <- matrix(c(s01.11.YY, s01.12.YS, s01.21.SY, s01.22.SS), nrow = 2, ncol = 2, byrow = TRUE)
  S.mat <- (1 / m.count) * S10 + (1 / n.count) * S01
  sd.y <- sqrt(S.mat[1, 1])
  sd.s <- sqrt(S.mat[2, 2])
  L <- t(as.matrix(c(1, -1)))

  sd.est <- (L %*% ((1 / m.count) * S10 + (1 / n.count) * S01) %*% t(L))^(1 / 2)

  return(list("u.y" = as.numeric(u.y), "u.s" = as.numeric(u.s), "delta.estimate" = as.numeric(delta.estimate), "sd.u.y" = sd.y, "sd.u.s" = sd.s, "sd.delta" = as.numeric(sd.est)))
}

test.surrogate <- function(full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, szero = NULL, epsilon = NULL, power.want.s = 0.7, u.y.hyp = NULL) {
  if (!is.null(full.data)) {
    yone <- full.data[full.data[, 3] == 1, 1]
    yzero <- full.data[full.data[, 3] == 0, 1]
    sone <- full.data[full.data[, 3] == 1, 2]
    szero <- full.data[full.data[, 3] == 0, 2]
  }
  dd <- delta.calculate(yone = yone, yzero = yzero, sone = sone, szero = szero)
  z.alpha <- qnorm(0.95)
  ci.delta <- c(-1, dd$delta.estimate + z.alpha * dd$sd.delta)

  if (is.null(epsilon)) {
    n1 <- length(yone)
    n0 <- length(yzero)
    sd.null <- sqrt((n1 + n0 + 1) / (12 * n1 * n0))
    z.alpha.2 <- qnorm(0.975)
    u.s.power <- 1 / 2 - (qnorm(1 - power.want.s) - z.alpha.2) * (sd.null)
    if (is.null(u.y.hyp)) {
      epsilon <- dd$u.y - u.s.power
    } else {
      epsilon <- u.y.hyp - u.s.power
    }
  }
  if (ci.delta[2] < epsilon) {
    is.surrogate <- TRUE
  }
  if (ci.delta[2] >= epsilon) {
    is.surrogate <- FALSE
  }

  return(list("u.y" = as.numeric(dd$u.y), "u.s" = as.numeric(dd$u.s), "delta.estimate" = as.numeric(dd$delta.estimate), "sd.u.y" = dd$sd.u.y, "sd.u.s" = dd$sd.u.s, "sd.delta" = as.numeric(dd$sd.delta), "ci.delta" = ci.delta, "epsilon.used" = epsilon, "is.surrogate" = is.surrogate))
}

est.power <- function(n.total, rho = 0.80, u.y.alt, delta.alt, power.want.s = 0.7) {
  delta.se.est <- sqrt(2 * (1 - rho) / (3 * n.total))
  sd.null <- sqrt((n.total + 1) / (3 * n.total^2))
  z.alpha.2 <- qnorm(0.975)
  u.s.power <- 1 / 2 - (qnorm(1 - power.want.s) - z.alpha.2) * (sd.null)
  z.alpha <- qnorm(0.95)
  est.power <- pnorm((1 / delta.se.est) * (u.y.alt - u.s.power - delta.alt) - z.alpha)
  return(est.power)
}

# a is one number, b is vector
var.wil <- function(a, b, flip = FALSE) {
  if (!flip) {
    return(mean(sapply(b, my.wilcox, x = a)))
  }
  if (flip) {
    return(mean(sapply(b, my.wilcox, y = a)))
  }
}

# in delong paper x = a, y = b
my.wilcox <- function(x, y) {
  if (x > y) {
    return(1)
  }
  if (x == y) {
    return(1 / 2)
  }
  if (x < y) {
    return(0)
  }
}
