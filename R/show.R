setMethod("show", "NormalPredInterval", function(object){
  cat("\n  Parametric Prediction Interval \nto contain at least", object@k, "out of", object@m, "future observations\n\n")
  cat("Level:", object@level, "\n")
  cat("Quantile:", round(object@quantile,3), "\n\n")  
  cat("   [", round(object@interval[1],3), ";", round(object@interval[2],3), "]\n\n")
})

setMethod("show", "nparPredInterval", function(object){
  cat("\n  Non-parametric Prediction Interval \nto contain at least", object@k, "out of", object@m, "future observations\n\n")
  cat("Level:", object@level, "\n\n")
  cat("   [", round(object@interval[1],3), ";", round(object@interval[2],3), "]\n\n")
})

setMethod("show", "NormalPrecInterval", function(object){
  cat("\n  Parametric Prediction Interval \nto contain the mean of", object@m, "future observations\n\n")
  cat("Level:", object@level, "\n\n")
  cat("   [", round(object@interval[1],3), ";", round(object@interval[2],3), "]\n\n")
})

setMethod("show", "nparPrecInterval", function(object){
  cat("\n  Non-parametric Prediction Interval \nto contain the median of", object@m, "future observations\n\n")
  cat("Level:", object@level, "\n\n")
  cat("   [", round(object@interval[1],3), ";", round(object@interval[2],3), "]\n\n")
})

setMethod("show", "PIlm", function(object){
  cat("\n  Prediction Interval \nto contain at least", object@k, "out of", nrow(object@newdata), "future observations\n\n")
  cat("Level:", object@level, "\n\n")
  print(object@prediction)
})