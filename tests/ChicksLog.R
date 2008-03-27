
library(richards)

library(nlme)

data(ChickWeight)

ChickWeight <- cbind(ChickWeight, ExpTime = exp(ChickWeight[,"Time"]))

groupedChicks <- groupedData(weight ~ ExpTime | Diet / Chick, 
                             inner = ~ PR, data = ChickWeight)

listChicks <- split(groupedChicks, list(groupedChicks[,"Diet"],
                                        groupedChicks[,"Chick"]),
                    drop = TRUE)

fit.RichardsLog.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSrichardsLog(Time, a, d, xmid, scal, ny),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.FplLog.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSnyFixedLog(1, Time, a, d, xmid, scal),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.FplLog.D <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSnyFixedLog(1, Time, a, d, xmid, scal),
          control = nls.control(warnOnly = TRUE), data = X))
 }

StartLog <- list(a =    25, d =  200, xmid =     1, scal =  0.1, ny =  1)
LowerLog <- list(a = -2000, d =    5, xmid =     0, scal =    0, ny =  0.01)
UpperLog <- list(a =    50, d = 1500, xmid =    20, scal =  100, ny = 128)

myApplyFUN2XLog <- 
  function(X, FUN, start = StartLog, lower = LowerLog, upper = UpperLog,
              Report = print, ...) {
    if (dim(X)[1] > 0) {
      result <- FUN(X, start, lower, upper)
      if (is.null(result) | class(result) == "try-error")
        Report(X)
      return(result)
    }
  }

fits.RichardsLog.BP <- fitList(listChicks, FUN = fit.RichardsLog.BP,
                               applyFUN2X = myApplyFUN2XLog, n = 5)

# x1.1 <- cbind(  fits.RichardsLog.BP[[3]],              fits.Richards.BP[[3]])
# x1.2 <- cbind(  fits.RichardsLog.BP[[3]][,"xmid"], log(fits.Richards.BP[[3]][,"x50"]))
# x1.3 <- cbind(1/fits.RichardsLog.BP[[3]][,"scal"],    (fits.Richards.BP[[3]][,"b"]))
# x1.4 <- cbind(  fits.RichardsLog.BP[[3]][,"scal"],  1/(fits.Richards.BP[[3]][,"b"]))

fits.FplLog.BP      <- fitList(listChicks, FUN = fit.FplLog.BP,
                               applyFUN2X = myApplyFUN2XLog, n = 4)

# x2.1 <- cbind(  fits.FplLog.BP[[3]],              fits.Fpl.BP[[3]])
# x2.2 <- cbind(  fits.FplLog.BP[[3]][,"xmid"], log(fits.Fpl.BP[[3]][,"x50"]))
# x2.3 <- cbind(1/fits.FplLog.BP[[3]][,"scal"],    (fits.Fpl.BP[[3]][,"b"]))
# x2.4 <- cbind(  fits.FplLog.BP[[3]][,"scal"],  1/(fits.Fpl.BP[[3]][,"b"]))

fits.FplLog.D       <- fitList(listChicks, FUN = fit.FplLog.D,
                               applyFUN2X = myApplyFUN2XLog, n = 4)

# x3.1 <- cbind(  fits.FplLog.D[[3]],              fits.Fpl.D[[3]])
# x3.2 <- cbind(  fits.FplLog.D[[3]][,"xmid"], log(fits.Fpl.D[[3]][,"x50"]))
# x3.3 <- cbind(1/fits.FplLog.D[[3]][,"scal"],    (fits.Fpl.D[[3]][,"b"]))
# x3.4 <- cbind(  fits.FplLog.D[[3]][,"scal"],  1/(fits.Fpl.D[[3]][,"b"]))

