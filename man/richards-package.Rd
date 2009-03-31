\name{richards-package}
\alias{richards-package}
%\alias{richards}
\docType{package}
\title{
  The package provides selfStarts and functions for Richards curves
}
\description{
 \link{selfStart} models and functions for Richards curves:

 \code{SSrichardsSimple <- selfStart( ~ d + (a - d) / (1 + (2^ny-1) * (input / x50)^b)^(1/ny), function (mCall, data, LHS) initialRichards(mCall, data, LHS), c("a", "d", "b", "x50", "ny"))}

 Why a package for this one line?

 The function for the initial values is a call to the \link{selfStart} model
 \link{SSfpl} of the four parameter logistic function.

 \link{selfStart} models for Richards curves with specific fixed values
 of the shape parameter \code{ny} are also given.
 For these a Richards function described by \code{ny}
 and the same asymptotes, location and slope parameters as
 for a fitted four parameter logistic function
 is found and used for the initial value.

 If \code{ny} is negative then \code{1 + (2^ny-1)*(input/x50)^b} can
 also be negative, and then this value cannot be raised to \code{1/ny}
 (unless \code{1/ny} is an integer).
 'Modified' \link{selfStart} models are given to handle this situation.

}
\details{
\tabular{ll}{
Package: \tab richards   \cr
Type:    \tab Package    \cr
Version: \tab 0.4.5      \cr
Date:    \tab 2008-09-01 \cr
License: \tab GPL Version 2 or later \cr
}
}
\author{
Jens Henrik Badsberg

Maintainer: Jens Henrik Badsberg <september2008richards@badsberg.eu>
}
\references{
Richards, F.J. 1959: A flexible growth function for empirical use.
J. Exp. Bot. 10: 290-300.

\code{http://en.wikipedia.org/wiki/Generalised_logistic_curve.}
}
\keyword{ package }
