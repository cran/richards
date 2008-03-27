
library(richards)

library(nlme)

demo(sourceRichards)


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


StartLog <- list(a =    25, d =  200, xmid =     1, scal =  0.1, ny =  1)
LowerLog <- list(a = -2000, d =    5, xmid =     0, scal =    0, ny =  0.01)
UpperLog <- list(a =    50, d = 1500, xmid =    20, scal =  100, ny = 128)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

fit.RichardsE.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSrichardsE(ExpTime, a, d, b, x50, ny),
                   ###########
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.RichardsE.BP <- fitList(listChicks, FUN = fit.RichardsE.BP,
            start = list(a =    25, d =  200, b =  10, e =    10, ny =  1),
            lower = list(a = -2000, d =    5, b =   0, e =     1, ny =  0.01),
            upper = list(a =    50, d = 1500, b =  50, e = 10^20, ny = 128),
                             applyFUN2X = myApplyFUN2X, n = 5)

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  #

fit.RichardsE.Raw.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ sSrichardsE(ExpTime, a, d, b, x50, ny),
                   ###########
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.RichardsE.Raw.BP <- fitList(listChicks, FUN = fit.RichardsE.Raw.BP,
            start = list(a =    25, d =  200, b =  10, e =    10, ny =  1),
            lower = list(a = -2000, d =    5, b =   0, e =     1, ny =  0.01),
            upper = list(a =    50, d = 1500, b =  50, e = 10^20, ny = 128),
                             applyFUN2X = myApplyFUN2X, n = 5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

fit.Richards.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSrichards(ExpTime, a, d, b, x50, ny),
                   ##########
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Richards.BP <- fitList(listChicks, FUN = fit.Richards.BP,
                            applyFUN2X = myApplyFUN2X, n = 5)

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  #

fit.Richards.Raw.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ sSrichards(ExpTime, a, d, b, x50, ny),
                   ##########
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.Richards.Raw.BP <- fitList(listChicks, FUN = fit.Richards.Raw.BP,
                            applyFUN2X = myApplyFUN2X, n = 5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

fit.RichardsLog.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSrichardsLog(Time, a, d, xmid, scal, ny),
                   #############
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.RichardsLog.BP <- fitList(listChicks, FUN = fit.RichardsLog.BP,
                               applyFUN2X = myApplyFUN2XLog, n = 5)

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  #

fit.RichardsLog.Raw.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ sSrichardsLog(Time, a, d, xmid, scal, ny),
                   #############
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.RichardsLog.Raw.BP <- fitList(listChicks, FUN = fit.RichardsLog.Raw.BP,
                               applyFUN2X = myApplyFUN2XLog, n = 5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #

fit.RichardsLogBM.C.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSrichardsLogBM(Time, a, d, xmid, b, m), # start = start, 
                   ###############
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.RichardsLogBM.C.BP <- fitList(listChicks, FUN = fit.RichardsLogBM.C.BP, n = 5,
            start = list(a =    25, d =  200, xmid =   1.00, b =   1.0, m =  1/1),
            lower = list(a = -2000, d =    5, xmid =   0.02, b =   0.0, m = -999),
            upper = list(a =    50, d = 1500, xmid =  20.00, b = 100.0, m =  999),
                              applyFUN2X = myApplyFUN2XLog)

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  #

fit.RichardsLogBM.Raw.C.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ sSrichardsLogBM(
                   ###############
          Time, a, d, xmid, b, m), # start = start, 
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fits.RichardsLogBM.Raw.C.BP <- fitList(listChicks, FUN = fit.RichardsLogBM.Raw.C.BP, n = 5,
            start = list(a =    25, d =  200, xmid =   1.00, b =   1.0, m =  1/1),
            lower = list(a = -2000, d =    5, xmid =   0.02, b =   0.0, m = -999),
            upper = list(a =    50, d = 1500, xmid =  20.00, b = 100.0, m =  999),
                              applyFUN2X = myApplyFUN2XLog)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  #
