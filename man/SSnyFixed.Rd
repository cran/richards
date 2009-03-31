\name{SSnyFixed}
\alias{SSnyFixed}
\alias{SSnyFixedLog}
\title{ Richards Function Model with the shape parameter 'ny' fixed }
\description{
This selfStart model evaluates a Richards function and its gradient. 
It has an initial attribute that will evaluate initial estimates
of the parameters \code{a}, \code{d}, \code{b},
and \code{b} for a given set of data for \code{ny} fixed. 

Here we fix the shape parameter \code{ny} to a specific value.
}
\usage{
           SSnyFixed(ny, input, a, d, b, x50)
        SSnyFixedLog(ny, input, a, d, xmid, scal)
}
\arguments{
    \item{ny}{ a numeric value. This has to be a 'number'!
               Sets the fixed shape parameter. }
 \item{input}{ a numeric vector of values at which to evaluate the model. }
     \item{a}{ a numeric parameter representing the horizontal asymptote on the
               left side (very small values of input) for \code{b} positive,
               else the horizontal asymptote on the right side. }
     \item{d}{ a numeric parameter representing the horizontal asymptote on the
               right side (very large values of input) for \code{b} positive,
               else the horizontal asymptote on the left side. }
     \item{b}{ a numeric scale parameter on the input axis,
               the 'growth rate', the reciprocal of the scale parameter
               \code{scal} for the four point logistic curve by \code{SSfpl}.
     (The 'growth' parameter \code{b} should be negative, and \code{a}
     is thus the right side (larger) asymptote, for Bertalenffy growth models.)
               }
   \item{x50}{ a numeric parameter representing the input value at
               the center of the curve:
               The value of \code{SSnyFixed} will be midway between
               \code{a} and \code{d} at \code{x50}. }
  \item{xmid}{ a numeric parameter representing the input value at
               the center of the curve:
               The value of \code{SSnyFixedLog} will be
               midway between \code{a} and \code{d} at \code{xmid}. }
  \item{scal}{ a numeric scale parameter on the input axis, the 'growth rate'. }
}
\details{
It works!
But this is probably not by the intended design of \code{selfStart}.

See \link{richards} about details of the Richard function.
Observe that the \code{selfStarts} models \code{SSnyfixed}
and \code{SSnyfixedLog} do not handle the Gompertz function.
See also \link{richards} about the modifications for
constant asymptote for \code{1 + (2^ny-1) * (input/x50)^b} or
\code{1 + (2^ny-1) * exp((input-xmid)/scal)} negative.
}
\value{
a numeric vector of the same length as \code{input}. 
It is the value of the expression
\code{d + (a - d) / (1 + (2^ny-1) * (input/x50)^b )^(1/ny)}, and
\code{d + (a - d) / (1 + (2^ny-1) * exp((input-xmid)/scal))^(1/ny)} for 
\code{SSnyFixed} and \code{SSnyFixedLog} respectively.
The value is \code{d} if \code{1 + (2^ny-1) * (input/x50)^b} or
\code{1 + (2^ny-1) * exp((input-xmid)/scal)} respectively is
less than \code{0}, else the above value.
If all of the arguments \code{a}, \code{d}, \code{b} (\code{scal}),
and \code{x50} (\code{xmid}) are names of objects,
the gradient matrix with respect to these names
is attached as an attribute named gradient.
}
\author{ Jens Henrik Badsberg }
\keyword{ package }
\examples{

# A fit of a richards function close to the Gompertz curve,
# but with negative 'ny':

fit.nyFixedLog.M0.01 <- function(X, lower, upper)
 if (dim(X)[1] > 0) {
  # Note: 'ny' as value, not the argument!
  result <- try(nls(SIGNAL ~ SSnyFixedLog(ny = -0.01, log(1 / CONC),
                                          a, d, xmid, scal),
                    lower = lower, upper = upper, data = X,
                    control = nls.control(warnOnly = TRUE), algorithm = "port"))
  result
 }

# If the fit gets close to the Gompertz curve,
# try to fit a Richards curve with negative 'ny'
# using the above function to get a starting value:

fit.RichardsLogBG.Double <- 
 function(X, start,
          # ArgLower and ArgUpper is assigned values when called from 'fitList':
          ArgLower = list(a = -0.5, d =  0.5, xmid = -6, b = 0.1, g = 0.01),
          ArgUpper = list(a =  0.5, d = 15.0, xmid =  6, b = 5.0, g = gmax),
          gmax = 100)
 if (dim(X)[1] > 0) {
  result <- NULL
  result <- try(nls(SIGNAL ~ SSrichardsLogBG(log(1/CONC), a, d, xmid, b, g),
                    lower = ArgLower, upper = ArgUpper, data = X,
                    control = nls.control(warnOnly = TRUE), algorithm = "port"))
  if (class(result) != "try-error") {
    parms <- summary(result)$parameters[, "Estimate"]; # print(c(1, parms))
    if (parms["g"] == gmax) {
      SlogLower <- list(a = -0.5, d =  0.5, scal = 0.2, xmid = -4)
      SlogUpper <- list(a =  0.5, d = 15.0, scal = 100, xmid =  4)
      result <- NULL
      result <- fit.nyFixedLog.M0.01(X, SlogLower, SlogUpper); # print(result)
      if (class(result) != "try-error") {
        parms <- summary(result)$parameters[, "Estimate"]; # print(c(2, parms))
        if (TRUE) {
          BlogLower <- list(a = -0.5, d =  0.5, xmid = -5, b = 0.1)
          BlogUpper <- list(a =  0.5, d = 15.0, xmid =  5, b = 5.0)
          result <- NULL
          result <- try(nls(SIGNAL ~ SSrichardsLogBG(log(1/CONC),
                                                     a, d, xmid, b, g),
                            start = list(a = parms["a"],
                                         d = parms["d"],
                                         b = 1 / parms["scal"],
                                         xmid = parms["xmid"],
                                         g = - gmax),
                            lower = append(BlogLower, list(g = -200)),
                            upper = append(BlogUpper, list(g = -0.4)),
                            control = nls.control(warnOnly = TRUE),
                            algorithm = "port", data = X)); # print(result)
        }
      }
    }
  }
  result
 }

# The data of the plates can not be supplied.
# Thus you can not run:

# fits.RichardsLogBM.Double.C <- fitList(listPlates,
#   FUN = fit.RichardsLogBM.Double, n = 5,
#   start = list(a =  0.0, d =  2.0, xmid =  0, b = 1.0, m = 0.5),
#   lower = list(a = -0.5, d =  0.5, xmid = -7, b = 0.1, m = 0.01),
#   upper = list(a =  0.5, d = 15.0, xmid =  7, b = 5.0, m = 100),
#                                   applyFUN2X = myApplyFUN2X)

}
