predint <- function(x, k, m, level=0.95, alternative="two.sided", quantile=NULL, absError=0.001, interval=c(0, 100)){
  if (!is.numeric(x)) stop("x is not a numeric vector!")
  k <- as.integer(k)
  m <- as.integer(m)
  if (k > m | k < 1) stop("k needs to be between 0 < k <= m")
  if (level < 0 | level > 1) stop("level has to be between 0 and 1")
  if (!alternative %in% c("two.sided", "less", "greater")) stop("Alternative has to be one of 'two.sided', 'less', or 'greater'")
  n <- length(x)
  if (is.null(quantile)){
    quant <- if (alternative == "two.sided") PIcritval(k, m, n, level, absError=absError, interval=interval) else PIonesided(k, m, n, level, absError=absError, interval=interval)
  } else quant <- quantile
  est <- mean(x)
  std <- sd(x)
  qterm <- quant * std * sqrt(1 + (1/n))
  lower <- if (alternative == "less") -Inf else est - qterm
  upper <- if (alternative == "greater") Inf else est + qterm
  new(Class="NormalPredInterval", quantile=quant, m=m, k=k, interval=c(lower, upper), sample=x, level=level, alternative=alternative)
}



nparpredint <- function(x, k, m, level=0.95, alternative="two.sided"){
  if (!is.numeric(x)) stop("x is not a numeric vector!")
  k <- as.integer(k)
  m <- as.integer(m)
  if (k > m | k < 1) stop("k needs to be between 0 < k <= m")
  if (level < 0 | level > 1) stop("level has to be between 0 and 1")
  if (!alternative %in% c("two.sided", "less", "greater")) stop("Alternative has to be one of 'two.sided', 'less', or 'greater'")  
  n <- length(x)
  dfl <- sapply(1:n, function(r) distfreetollim(r=r, m=m, n=n, K=level))
  nu <- sum(dfl >= k)
  if (nu == 0){
    lower <- -Inf
    upper <- Inf
    warning("The prediction limits are set to min(x)/max(x). The interval may be inaccurate!")
  } else {
    if (alternative == "two.sided"){
      v <- c(floor(nu/2), ceiling(nu/2))
      v1 <- v2 <- v[v > 0]
      lower <- mean(x[order(x)][v1])
      upper <- mean(x[order(x)][n-v2+1])
    }
    if (alternative == "lower"){
      lower <- x[order(x)][nu]
      upper <- Inf
    }
    if (alternative == "upper"){
      lower <- -Inf
      upper <- x[order(x)][n-nu+1]
    }
  }
  new(Class="nparPredInterval", m=m, k=k, interval=c(lower, upper), sample=x, level=level, alternative=alternative)
}




