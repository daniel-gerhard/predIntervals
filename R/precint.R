precint <- function(x, m, level=0.95, alternative="two.sided"){
  if (!is.numeric(x)) stop("x is not a numeric vector!")
  m <- as.integer(m)
  if (level < 0 | level > 1) stop("level has to be between 0 and 1")
  if (!alternative %in% c("two.sided", "less", "greater")) stop("Alternative has to be one of 'two.sided', 'less', or 'greater'")
  n <- length(x)
  quant <- if (alternative == "two.sided") qt(1-(1-level)/2, df=n-1) else qt(level, df=n-1)
  est <- mean(x)
  std <- sd(x)
  qterm <- quant * std * sqrt(1/m + (1/n))
  lower <- if (alternative == "less") -Inf else est - qterm
  upper <- if (alternative == "greater") Inf else est + qterm
  new(Class="NormalPrecInterval", m=m, interval=c(lower, upper), sample=x, level=level, alternative=alternative)
}


nparprecint <- function(x, m, level=0.95, alternative="two.sided"){
  if (!is.numeric(x)) stop("x is not a numeric vector!")
  m <- as.integer(m)
  if (level < 0 | level > 1) stop("level has to be between 0 and 1")
  if (!alternative %in% c("two.sided", "less", "greater")) stop("Alternative has to be one of 'two.sided', 'less', or 'greater'")
  n <- length(x)
  conf.level <- if (alternative == "two.sided") level else 1-2*(1-level)
  j <- median(1:m)
  w <- 1:n
  ap <- (choose(j+w-1,w) * choose(m+n-j-w,n-w)) / choose(m+n,n)
  a <- 0
  b <- n-a
  while (sum(ap[a:b]) > conf.level){
    a <- a+1
    b <- n-a
  }
  a <- a-1
  b <- n-a+1
  lower <- if (alternative == "less") -Inf else x[order(x)][a]
  upper <- if (alternative == "greater") Inf else x[order(x)][b]
  new(Class="nparPrecInterval", m=m, interval=c(lower, upper), sample=x, level=level, alternative=alternative)
}

