setClass("PInt", representation(interval = "numeric",
                                sample = "numeric",
                                level = "numeric",
                                alternative = "character"))


setClass("PredInterval", representation(m = "numeric",
                                        k = "numeric"),
         contains="PInt")
setClass("NormalPredInterval", representation(quantile = "numeric"),
         contains="PredInterval")
setClass("nparPredInterval", contains="PredInterval")


setClass("PrecInterval", representation(m = "numeric"),
         contains="PInt")
setClass("NormalPrecInterval", contains="PrecInterval")
setClass("nparPrecInterval", contains="PrecInterval")

setClass("PIlm", representation(prediction = "data.frame", 
                                newdata = "data.frame",
                                model = "lm",
                                k = "numeric",
                                level = "numeric",
                                alternative = "character"))