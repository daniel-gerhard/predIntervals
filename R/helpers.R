PIcritval <- function(k,m,n,alpha,nu, absError=0.001, interval=c(0, 100)){
  p <- 1/(n + 1)
  intfunc <- function(x, j, u, nu, n, m, p){
    t <- x[1] / (1-x[1]^2)
    s <- x[2] / (1-x[2])
    a <- (-u * s + (p)^(0.5) * t)/((1 - p)^(0.5))
    b <- (u * s + (p)^(0.5) * t)/((1 - p)^(0.5))
    w <- pnorm(b) - pnorm(a)
    r1 <- choose(m, j) * (w)^j * (1 - w)^(m - j) * dnorm(t) 
    fnu <- (2 * nu * s * dchisq(nu * s^2, nu)) 
    f1 <- (1 + x[1]^2) / ((1-x[1]^2)^2) 
    f2 <- 1 / ((1-x[2])^2) 
    out <- r1 * fnu * f1 * f2
    return(out)
  }
  adaptintfunc <- function(j, u) adaptIntegrate(intfunc, lowerLimit=c(-1, 0), upperLimit=c(1, 1), j=j, u=u, nu=nu, n=n, m=m, p=p, absError=absError)$integral
  helper <- function(u) {
    sapply(u, function(u) {
      (sum(sapply(k:m, function(j) adaptintfunc(j, u)))) - alpha
    })
  }
  ustar <- uniroot(helper, interval)$root
  return(ustar)  
}


##########################################

PIonesided <- function(k,m,n,alpha,nu, absError=0.001, interval=c(0, 100)){
  p <- 1/(n + 1)
  intfunc <- function(x, j, u, nu, n, m, p){
    t <- x[1] / (1-x[1]^2)
    s <- x[2] / (1-x[2])
    b <- (u * s + (p)^(0.5) * t)/((1 - p)^(0.5))
    w <- pnorm(b)
    r1 <- choose(m, j) * (w)^j * (1 - w)^(m - j) * dnorm(t) 
    fnu <- (2 * nu * s * dchisq(nu * s^2, nu)) 
    f1 <- (1 + x[1]^2) / ((1-x[1]^2)^2) 
    f2 <- 1 / ((1-x[2])^2) 
    out <- r1 * fnu * f1 * f2
    return(out)
  }
  adaptintfunc <- function(j, u) adaptIntegrate(intfunc, lowerLimit=c(-1, 0), upperLimit=c(1, 1), j=j, u=u, nu=nu, n=n, m=m, p=p, absError=absError)$integral
  helper <- function(u) {
    sapply(u, function(u) {
      (sum(sapply(k:m, function(j) adaptintfunc(j, u)))) - alpha
    })
  }
  ustar <- uniroot(helper, interval)$root
  return(ustar)  
}


##########################################

distfreetollim <- function(r, m, n, K){
  N0 <- 0:m
  PN0 <- (choose(N0 + n - r, N0) * choose(m - N0 + r - 1, m - N0)) / (choose(n + m, m))
  PN0i <- PN0[(m+1):1]
  sum(cumsum(PN0i) >= K)-1
}

