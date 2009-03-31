\name{llamodified}
% \alias{fitList}
% \alias{fitDoubleList}
\alias{lamodified}
\alias{llamodified}

\alias{amList-class}
\alias{amLlist-class}

\alias{coef,amList-method}
\alias{show,amList-method}
\alias{summary,amList-method}

\alias{coef,amLlist-method}
\alias{show,amLlist-method}
\alias{summary,amLlist-method}

\title{ Apply a Function (after modifying argument) over a List or List of lists }

\description{ 

\code{lamodified} returns a list of the same length as \code{list}, 
each element of which is the result of applying \code{applyFUN2X}
to the corresponding element of \code{list}. 

\code{llamodified} works similar, but on a list of lists.
}

\usage{
   lamodified(list, FUN,
              applyFUN2X = function(X, FUN, ...) FUN(X, ...), ...) 
  llamodified(list, FUN, 
              applyFUN2X = function(X, FUN, ...) FUN(X, ...), ...) 

%       fitList(list, FUN,
%               applyFUN2X = function(X, FUN, ...) FUN(X, ...),
%               n = 4, ...)
% fitDoubleList(list, FUN,
%               applyFUN2X = function(X, FUN, ...) FUN(X, ...),
%               n = 4, ...)

}

\arguments{
        \item{list}{ a list of datasets. }
         \item{FUN}{ a function estimating a model for each dataset. }
  \item{applyFUN2X}{ a function applied to each dataset,
                     retrieving starting value and fixing outliers,
                     and then calling the function \code{FUN} to
                     estimate the model for the modified dataset. }
%          \item{n}{ a numeric with the number of parameters of the model. }
       \item{\dots}{ additional arguments to \code{applyFUN2X}. }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "amlist")} and
                \code{signature(object = "amLlist")}: displays fits
                (if \code{FUN} returns an object with the method \code{coef}). }
    \item{coef}{\code{signature(object = "amlist")} and
                \code{signature(object = "amLlist")}: return parameters
                (if \code{FUN} returns an object with the method \code{coef}. }
    \item{summary}{\code{signature(object = "amlist")} and
                   \code{signature(object = "amLlist")}: return parameters with tests
                   (if \code{FUN} returns an object with the method \code{summary}. }
  }
}

\details{ 
a list of objects of class returned by (\code{FUN}) with
as many components as the number of item in the argument list.
Generic functions such as \code{coef}, \code{summary}, and \code{show}
have methods that can be applied to an \code{amList} or \code{amLlist} object.

Alternatives to \code{amList} are \link[base]{lapply},
\link[nlme]{lmList} and \link[nlme]{nlsList}.
In difference to \link[nlme]{lmList} and \link[nlme]{nlsList}
data is grouped before applying \code{amList}, since
the grouped data together with fits is necessary for backfitted values.

\code{llamodified} (\code{fitDoubleList}) is for 'similarity',
a list of 'experiments', each with a list of 'samples'.

% \code{fitList} and \code{fitDoubleList} are depricated.
% They will olne by retained until I have updated the examples, demos, and test-runs.

}

\author{ Jens Henrik Badsberg }

\keyword{ package }
\examples{

## Not run: 

library(richards)
library(nlme)

data(simulatedPlates)

groupedData <- groupedData(SIGNAL ~ CONC | PT / LAY / PLA, 
                           inner = ~ PR, data = simulatedPlates)

fullData <- split(groupedData, list(groupedData[,"PT"],
                                    groupedData[,"LAY"],
                                    groupedData[,"PLA"]),
                  drop = TRUE)

n <- length(fullData)
listData <- fullData[1:28]

fit.Richards.BP <- function(X, start, lower, upper)
 {
  try(nls(SIGNAL ~ SSrichards(CONC, a, d, b, x50, ny),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(maxiter = 1000, tol = 1e-04, 
                                warnOnly = TRUE), data = X))
 }

fit.Gompertz.BP <- function(X, start, lower, upper)
 {
  try(nls(SIGNAL ~ SSny0(CONC, a, d, b, x50),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.GompertzLog.D <- function(X, start, lower, upper)
 {
  try(nls(SIGNAL ~ SSny0Log(log(CONC), a, d, xmid, scal),
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.GompertzLogB.D <- function(X, start, lower, upper)
 {
  try(nls(SIGNAL ~ SSny0LogB(log(CONC), a, d, xmid, b),
          control = nls.control(warnOnly = TRUE), data = X))
 }

fit.Fpl.BP <- function(X, start, lower, upper)
 {
  try(nls(SIGNAL ~ SSny1(CONC, a, d, b, x50),
          lower = lower, upper = upper, algorithm = "port",
          control = nls.control(maxiter = 1000, tol = 1e-04, 
                                warnOnly = TRUE), data = X))
 }

fit.Fpl.D <- function(X, start, lower, upper)
 {
  try(nls(SIGNAL ~ SSny1(CONC, a, d, b, x50),
          control = nls.control(warnOnly = TRUE), data = X))
 }

Start <- list(a =  0.0, d =  3.0,  b = 1.0, x50 = 0.10,  ny =   0.1)
Lower <- list(a = -0.5, d =  0.02, b = 0.1, x50 = 0.001, ny =  -5)
Upper <- list(a =  0.5, d = 15.0,  b = 2.5, x50 = 5.00,  ny = 128)

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

fits.Richards.BP    <- lamodified(listData, FUN = fit.Richards.BP,
                                  applyFUN2X = myApplyFUN2X, n = 5)
coef(fits.Richards.BP)

fits.Gompertz.BP    <- lamodified(listData, FUN = fit.Gompertz.BP,
                                  applyFUN2X = myApplyFUN2X, n = 4)
coef(fits.Gompertz.BP)

fits.GompertzLog.D <- lamodified(listData, FUN = fit.GompertzLog.D,
                                 applyFUN2X = myApplyFUN2X, n = 4)
coef(fits.GompertzLog.D)

fits.GompertzLogB.D <- lamodified(listData, FUN = fit.GompertzLogB.D,
                                  applyFUN2X = myApplyFUN2X, n = 4)
coef(fits.GompertzLogB.D)

fits.Fpl.BP         <- lamodified(listData, FUN = fit.Fpl.BP,
                                  applyFUN2X = myApplyFUN2X, n = 4)
coef(fits.Fpl.BP)

fits.Fpl.D          <- lamodified(listData, FUN = fit.Fpl.D,
                                  applyFUN2X = myApplyFUN2X, n = 4)
coef(fits.Fpl.D)

X <- cbind(coef(fits.Richards.BP), coef(fits.Fpl.BP))

## End(Not run)

}
