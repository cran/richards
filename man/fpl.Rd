\name{fpl}
\alias{fpl}
\alias{fpl.deriv}
\alias{fpl.deriv.x}
\alias{fpl.inv}
\title{ Four parameter logistic function }
\description{ The four parameter logistic function, the inverse of the
 function, and the derivative of \code{fpl(x)} and \code{fpl(exp(x))}. }
\usage{
        fpl(x, a = 0.1, d = 2.4, e = 100, b = 1)
  fpl.deriv(x, a = 0.1, d = 2.4, e = 100, b = 1)
fpl.deriv.x(x, a = 0.1, d = 2.4, e = 100, b = 1)
    fpl.inv(y, a = 0.1, d = 2.4, e = 100, b = 1)
}
\arguments{
  \item{x}{ a numeric vector of values at which to evaluate the function. }
  \item{y}{ a numeric vector of values at which to evaluate the inverse function. }
  \item{a}{ a numeric parameter representing the horizontal asymptote on the
            left side (very small values of input) for \code{b} positive,
            else the horizontal asymptote on the right side. }
  \item{d}{ a numeric parameter representing the horizontal asymptote on the
            right side (very large values of input) for \code{b} positive,
            else the horizontal asymptote on the left side. }
  \item{b}{ a numeric scale parameter on the input axis, 'growth rate'.  }
  \item{e}{ the parameter representing the input value at the 
            inflection point of the curve. The value of \code{fpl}
            will be midway between \code{a} and \code{d} at \code{e}. }
}
\value{ A vector with function values. }
\author{ Jens Henrik Badsberg }
\keyword{ package }
