library(richards)

library(nlme)

data(ChickWeight)

ChickWeight <- cbind(ChickWeight, ExpTime = exp(ChickWeight[,"Time"]))

groupedChicks <- groupedData(weight ~ ExpTime | Diet / Chick, 
                             inner = ~ PR, data = ChickWeight)

listChicks <- split(groupedChicks, list(groupedChicks[,"Diet"],
                                        groupedChicks[,"Chick"]),
                    drop = TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

Start <- list(a =    25, d =  200, b =  10, x50 =    10, ny =  1)
Lower <- list(a = -2000, d =    5, b =   0, x50 =     1, ny =  0.01)
Upper <- list(a =    50, d = 1500, b =  50, x50 = 10^20, ny = 128)

myApplyFUN2X <- 
  function(X, FUN, start = Start, lower = Lower, upper = Upper,
              Report = print, ...) {
    if (dim(X)[1] > 0) {
      result <- FUN(X, start, lower, upper)
      # if (is.null(result) | class(result) == "try-error")
      #  Report(X)
      return(result)
    }
  }


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  #

demo(sourceNyExtras)

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  #

fit.Fpl.half.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ sSny.1o2(ExpTime, a, d, b, x50),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.Fpl.half.D <- function(X, start, lower, upper)
 {
  try(nls(weight ~ sSny.1o2(ExpTime, a, d, b, x50),
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Fpl.half.BP <- fitList(listChicks, FUN = fit.Fpl.half.BP,
                             applyFUN2X = myApplyFUN2X, n = 4)

fits.Fpl.half.D  <- fitList(listChicks, FUN = fit.Fpl.half.D,
                             applyFUN2X = myApplyFUN2X, n = 4)


fit.Fpl.third.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny.1o3(ExpTime, a, d, b, x50),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Fpl.third.BP <- fitList(listChicks, FUN = fit.Fpl.third.BP,
                             applyFUN2X = myApplyFUN2X, n = 4)

