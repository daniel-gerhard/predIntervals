lmpredint <- function(object, newdata, k, level = 0.95, alternative="two.sided", quantile=NULL, absError=0.001, interval=c(0, 100)){
  tt <- terms(object)
  if (!inherits(object, "lm")) warning("calling predict.lm(<fake-lm-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    stop("newdata is missing!")
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.pass, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }
  if (k > nrow(newdata) | k < 1) stop("k needs to be between 0 < k <= nrow(newdata)")
  n <- length(object$residuals)
  p <- object$rank
  p1 <- seq_len(p)
  piv <- if (p) stats:::qr.lm(object)$pivot[p1]
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) warning("prediction from a rank-deficient fit may be misleading")
  beta <- object$coefficients
  predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
  if (!is.null(offset)) predictor <- predictor + offset
  
  w <- stats:::weights.default(object)
  if (!is.null(w)) {
    weights <- w
    warning("assuming prediction variance inversely proportional to weights used for fitting\n")
  }
  
  weights <- 1
  if (!is.null(object$weights)) warning("Assuming constant prediction variance even though model fit is weighted\n")
  if (inherits(weights, "formula")) {
    if (length(weights) != 2L) stop("'weights' as formula should be one-sided")
    d <- newdata
    weights <- eval(weights[[2L]], d, environment(weights))
  }
  
  
  r <- object$residuals
  w <- object$weights
  rss <- sum(if (is.null(w)) r^2 else r^2 * w)
  df <- object$df.residual
  res.var <- rss/df
  
  if (p > 0) {
    XRinv <- X[, piv] %*% qr.solve(qr.R(stats:::qr.lm(object))[p1, p1])
    ip <- drop(XRinv^2 %*% rep(res.var, p))
  } else ip <- rep(0, n)
  
  if (is.null(quantile)){
    if (alternative == "two.sided"){
      tfrac <- -PIlmcritval(k,m=nrow(newdata),n,level, df, absError=absError, interval=interval)
    } else {
      tfrac <- -PIlmonesided(k,m=nrow(newdata),n,level, df, absError=absError, interval=interval)
    }
  } else {
    tfrac <- -quantile
  }
  
  pred.var = res.var/weights
  hwid <- tfrac * sqrt(ip + pred.var)
  predictor <- data.frame(cbind(predictor, predictor + hwid %o% c(1, -1)))
  colnames(predictor) <- c("fit", "lower", "upper")
  if (alternative == "less") predictor[,2] <- -Inf
  if (alternative == "greater") predictor[,3] <- Inf
  
  new(Class="PIlm", 
      prediction = predictor, 
      newdata = newdata,
      model = object,
      k = k,
      level = level,
      alternative = alternative)
}
