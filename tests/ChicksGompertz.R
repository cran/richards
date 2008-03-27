
library(richards)

library(nlme)

data(ChickWeight)

ChickWeight <- cbind(ChickWeight, ExpTime = exp(ChickWeight[,"Time"]))

groupedChicks <- groupedData(weight ~ ExpTime | Diet / Chick, 
                             inner = ~ PR, data = ChickWeight)

listChicks <- split(groupedChicks, list(groupedChicks[,"Diet"],
                                        groupedChicks[,"Chick"]),
                    drop = TRUE)

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Ny0

fit.Ny0.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny0(ExpTime, a, d, b, x50),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Ny0.BP      <- fitList(listChicks, FUN = fit.Ny0.BP,
                            applyFUN2X = myApplyFUN2X, n = 4)

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . #

fit.Ny0.D <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny0(ExpTime, a, d, b, x50),
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Ny0.D       <- fitList(listChicks, FUN = fit.Ny0.D,
                            applyFUN2X = myApplyFUN2X, n = 4)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# SSny0LogB

fit.Ny0LogB.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny0LogB(Time, a, d, xmid, b),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Ny0LogB.BP      <- fitList(listChicks, FUN = fit.Ny0LogB.BP,
                                applyFUN2X = myApplyFUN2XLog, n = 4)

cbind(  fits.Ny0LogB.BP[[3]],              fits.Ny0.BP[[3]])
cbind(  fits.Ny0LogB.BP[[3]][,"xmid"], log(fits.Ny0.BP[[3]][,"x50"]))
cbind(  fits.Ny0LogB.BP[[3]][,"b"],       (fits.Ny0.BP[[3]][,"b"]))

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . #

fit.Ny0LogB.D <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny0LogB(Time, a, d, xmid, b),
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Ny0LogB.D       <- fitList(listChicks, FUN = fit.Ny0LogB.D,
                                applyFUN2X = myApplyFUN2XLog, n = 4)

cbind(  fits.Ny0LogB.D[[3]],              fits.Ny0.D[[3]])
cbind(  fits.Ny0LogB.D[[3]][,"xmid"], log(fits.Ny0.D[[3]][,"x50"]))
cbind(  fits.Ny0LogB.D[[3]][,"b"],       (fits.Ny0.D[[3]][,"b"]))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# SSny0Log

fit.Ny0Log.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny0Log(Time, a, d, xmid, scal),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Ny0Log.BP      <- fitList(listChicks, FUN = fit.Ny0Log.BP,
                               applyFUN2X = myApplyFUN2XLog, n = 4)

cbind(  fits.Ny0Log.BP[[3]],              fits.Ny0.BP[[3]])
cbind(  fits.Ny0Log.BP[[3]][,"xmid"], log(fits.Ny0.BP[[3]][,"x50"]))
cbind(1/fits.Ny0Log.BP[[3]][,"scal"],    (fits.Ny0.BP[[3]][,"b"]))

cbind(fits.Ny0LogB.BP[[3]],            fits.Ny0Log.BP[[3]])
cbind(fits.Ny0LogB.BP[[3]][,"xmid"],   fits.Ny0Log.BP[[3]][,"xmid"])
cbind(fits.Ny0LogB.BP[[3]][,"b"],    1/fits.Ny0Log.BP[[3]][,"scal"])

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . #

fit.Ny0Log.D <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny0Log(Time, a, d, xmid, scal),
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Ny0Log.D       <- fitList(listChicks, FUN = fit.Ny0Log.D,
                               applyFUN2X = myApplyFUN2XLog, n = 4)

cbind(  fits.Ny0Log.D[[3]],              fits.Ny0.D[[3]])
cbind(  fits.Ny0Log.D[[3]][,"xmid"], log(fits.Ny0.D[[3]][,"x50"]))
cbind(1/fits.Ny0Log.D[[3]][,"scal"],    (fits.Ny0.D[[3]][,"b"]))

cbind(fits.Ny0LogB.D[[3]],            fits.Ny0Log.D[[3]])
cbind(fits.Ny0LogB.D[[3]][,"xmid"],   fits.Ny0Log.D[[3]][,"xmid"])
cbind(fits.Ny0LogB.D[[3]][,"b"],    1/fits.Ny0Log.D[[3]][,"scal"])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Gompertz

fit.gompertz.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSgompertz(Time, d, xmid, scal),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.gompertz.BP      <- fitList(listChicks, FUN = fit.gompertz.BP,
                                 applyFUN2X = myApplyFUN2XLog, n = 3)

x5.1 <- cbind(  fits.Ny0Log.D[[3]],   fits.Ny0.D[[3]])
x5.2 <- cbind(  fits.Ny0Log.BP[[3]],  fits.Ny0.BP[[3]])
x5.3 <- cbind(  fits.Ny0Log.BP[[3]],  fits.Ny0.BP[[3]], fits.gompertz.BP[[3]])
x5.4 <- cbind(  fits.Ny0Log.BP[[3]],  fits.Ny0.BP[[3]], fits.gompertz.BP[[3]])[,c(2,6,9)]

# format(x5.4, scientific = FALSE)

#
