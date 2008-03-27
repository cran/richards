\name{fitList}
\alias{fitList}
\alias{fitDoubleList}
\title{ Fit Richard curves to a list of data sets }
\description{ Fit Richard curves to a list of data sets. }
\usage{
      fitList(list, FUN,
              applyFUN2X = function(X, FUN, ...) FUN(X, ...),
              n = 4, ...)
fitDoubleList(list, FUN,
              applyFUN2X = function(X, FUN, ...) FUN(X, ...),
              n = 4, ...)
}
\arguments{
        \item{list}{ a list of datasets. }
         \item{FUN}{ a function estimating a model for each dataset. }
  \item{applyFUN2X}{ a function applied to each dataset,
                     retrieving starting value and fixing outliers,
                     and then calling the function \code{FUN} to
                     estimate the model for the modified dataset. }
           \item{n}{ a numeric with the number of parameters of the model. }
       \item{\dots}{ additional arguments to \code{applyFUN2X}. }
}
\details{ 'fitDoubleList' is for 'similarity', a list of 'experiments',
 each with a list of 'samples'. }
\author{ Jens Henrik Badsberg }
\keyword{ package }
\examples{

## Not run: 

library(nlme)

data(ChickWeight)

ChickWeight <- cbind(ChickWeight, ExpTime = exp(ChickWeight[,"Time"]))

groupedChicks <- groupedData(weight ~ ExpTime | Diet / Chick, 
                             inner = ~ PR, data = ChickWeight)

fullChicks <- split(groupedChicks, list(groupedChicks[,"Diet"],
                                        groupedChicks[,"Chick"]),
                    drop = TRUE)

n <- length(fullChicks)
listChicks <- fullChicks[(n-14):n]

fit.Richards.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSrichards(ExpTime, a, d, b, x50, ny),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.Fpl.BP <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny1(ExpTime, a, d, b, x50),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.Fpl.D <- function(X, start, lower, upper)
 {
  try(nls(weight ~ SSny1(ExpTime, a, d, b, x50),
          control = nls.control(warnOnly = TRUE), data = X))
 }

Start <- list(a =    25, d =  200, b = 10, x50 =    10, ny =  1)
Lower <- list(a = -2000, d =    5, b =  0, x50 =     1, ny =  0.01)
Upper <- list(a =    50, d = 1500, b = 50, x50 = 10^20, ny = 128)

myApplyFUN2X <- 
  function(X, FUN, start = Start, lower = Lower, upper = Upper,
              Report = print, ...) {
    if (dim(X)[1] > 0) {
      result <- FUN(X, start, lower, upper)
      if (is.null(result) | class(result) == "try-error")
        Report(X)
      return(result)
    }
  }

fits.Richards.BP <- fitList(listChicks, FUN = fit.Richards.BP,
                            applyFUN2X = myApplyFUN2X, n = 5)
fits.Fpl.BP      <- fitList(listChicks, FUN = fit.Fpl.BP,
                            applyFUN2X = myApplyFUN2X, n = 4)
fits.Fpl.D       <- fitList(listChicks, FUN = fit.Fpl.D,
                            applyFUN2X = myApplyFUN2X, n = 4)

cbind(fits.Richards.BP$estimates, fits.Fpl.BP$estimates)

## End(Not run)

}
