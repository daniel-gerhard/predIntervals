PIcritval <- function(k,m,n,alpha){
  nu <- n - 1
  p <- 1/(n + 1)
  fnu <- function(s) {
    2 * nu * s * dchisq(nu * s^2, n-1)
  }
  g <- function(s, j, u) {
    integrate(function(t) {
      sapply(t, function(t) {
        a <- (-u * s + (p)^(0.5) * t)/((1 - p)^(0.5))
        b <- (u * s + (p)^(0.5) * t)/((1 - p)^(0.5))
        w <- pnorm(b) - pnorm(a)
        choose(m, j) * (w)^j * (1 - w)^(m - j) * dnorm(t)
      })
    }, -Inf, Inf)$value * fnu(s)
  }
  PB <- function(j, u) {
    integrate(function(s) {
      sapply(s, function(s) {
        g(s, j, u)
      })
    }, 0, Inf)$value
  }
  h <- function(u) {
    sapply(u, function(u) {
      (sum(sapply(k:m, function(j) PB(j, u)))) - alpha
    })
  }
  ustar <- uniroot(h, c(0, 100))$root
  sqrt((n+1)/n)*ustar
}

##########################################

PIonesided <- function(k,m,n,alpha){
  nu <- n - 1
  p <- 1/(n + 1)
  fnu <- function(s) {
    2 * nu * s * dchisq(nu * s^2, n-1)
  }
  g <- function(s, j, u) {
    integrate(function(t) {
      sapply(t, function(t) {
        b <- (u * s + (p)^(0.5) * t)/((1 - p)^(0.5))
        w <- pnorm(b)
        choose(m, j) * (w)^j * (1 - w)^(m - j) * dnorm(t)
      })
    }, -Inf, Inf)$value * fnu(s)
  }
  PB <- function(j, u) {
    integrate(function(s) {
      sapply(s, function(s) {
        g(s, j, u)
      })
    }, 0, Inf)$value
  }
  h <- function(u) {
    sapply(u, function(u) {
      (sum(sapply(k:m, function(j) PB(j, u)))) - alpha
    })
  }
  ustar <- uniroot(h, c(0, 100))$root
  sqrt((n+1)/n)*ustar
}

##########################################

distfreetollim <- function(r, m, n, K){
  N0 <- 0:m
  PN0 <- (choose(N0 + n - r, N0) * choose(m - N0 + r - 1, m - N0)) / (choose(n + m, m))
  PN0i <- PN0[(m+1):1]
  sum(cumsum(PN0i) >= K)-1
}

