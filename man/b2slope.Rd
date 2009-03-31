\name{b2slope}
\alias{b2slope}
\alias{slope2b}
\alias{ny2f}
\alias{x2c}
\title{The slope and the parameter 'b' of Richards curves}
\description{
Finds the slope at \code{x50} of a Richards curve,
or the parameter \code{b} of a Richards curve with a given slope.
\code{x2c} is for finding the parameters \code{C} and \code{E}
from \code{x50}.
}
\usage{
b2slope(ny = X[, "ny"], b = X[, "b"], a = X[, "a"], d = X[, "d"], 
        X = NULL, fold = 2, delta = 1e-09)
slope2b(ny = X[, "ny"], slope = X[, "slope"], a = X[, "a"], d = X[, "d"], 
        X = NULL, fold = 2, delta = 1e-09)
ny2f(ny, delta = 1e-09)
x2c(x50 = X[, "x50"], ny = X[, "ny"], b = X[, "b"],
    X = NULL, fold = 2, x50f = 1, delta = 1e-09)
}
\arguments{
    \item{ny}{ the shape parameter: affects near which asymptote maximum
               'growth' occurs. }
     \item{b}{ a numeric scale parameter on the input axis,
               the 'growth rate', the reciprocal of the scale parameter
               \code{scal} for the four point logistic curve by SSfpl. 
     (The 'growth' parameter \code{b} should be negative, and \code{a}
     is thus the right side (larger) asymptote, for Bertalenffy growth models.)
              }
 \item{slope}{ a numeric parameter representing the slope at
               \code{x50} of the curve. }
     \item{a}{ a numeric parameter representing the horizontal asymptote on the
               left side (very small values of input) for \code{b} positive,
               else the horizontal asymptote on the right side. }
     \item{d}{ a numeric parameter representing the horizontal asymptote on the
               right side (very large values of input) for \code{b} positive,
               else the horizontal asymptote on the left side. }
   \item{x50}{ a numeric parameter representing the input value at
               the center of the curve:
               The value of \code{SSrichards} will be
               midway between \code{a} and \code{d} at \code{x50}. }
  \item{x50f}{ a numeric value, factor for \code{x50}. }
     \item{X}{ a numeric matrix with the parameters of the Richards curves. }
  \item{fold}{ a numeric with the base of the logarithm. }
 \item{delta}{ a numeric value.
               If \code{abs(ny) < delta} then the Gompertz curve is used. }
}
\details{
The slope (with respect to the logarithm with base \code{fold} of \code{x})
of the Richards curve is \code{b * (d-a) * f * log(fold)},
with \code{f} given by \code{ny2f}.
}
\value{A numeric vector}
\author{ Jens Henrik Badsberg }
\keyword{ package }
