\name{richards}
\alias{richards}
\alias{richards.deriv}
\alias{richards.deriv.x}
\alias{richards.inv}
\alias{solveB}
\alias{solveE}
\title{ Richards functions: Five parameter (5PL) logistic curves }
\description{ The Richards function, the inverse of
the function, and the derivative of \code{richards(x)}
(and the derivative of \code{richards(x)} with
respect to \code{log(x)} in the point \code{x},
that is derivative of \code{richards(exp(u))} with
respect to \code{u} in the point \code{x = exp(u)}). 

In other implementation of Richards functions and
four parameter logistic functions the argument is
\code{log(x)} of the here used argument \code{x}.

Here we have a positive concentration in mind,
where value depending on the concentration can
be described by Richards functions
in the logarithm of the concentration.
}
\usage{
        richards(x, a = 0.1, d = 2.4,
                 e = solveE(x50, b, ny), x50 = 100,
                 b = solveB(x = x50, a = a, d = d, ny = ny, 
                            x50 = x50, b4 = b4), b4 = 1,
                 ny = k - 1, k = 2)

  richards.deriv(x, a = 0.1, d = 2.4,
                 e = solveE(x50, b, ny), x50 = 100,
                 b = solveB(x = x50, a = a, d = d, ny = ny, 
                            x50 = x50, b4 = b4), b4 = 1,
                 ny = k - 1, k = 2)

richards.deriv.x(x, a = 0.1, d = 2.4,
                 e = solveE(x50, b, ny), x50 = 100,
                 b = solveB(x = x50, a = a, d = d, ny = ny, 
                            x50 = x50, b4 = b4), b4 = 1,
                 ny = k - 1, k = 2)

    richards.inv(y, a = 0.1, d = 2.4,
                 e = solveE(x50, b, ny), x50 = 100,
                 b = solveB(x = x50, a = a, d = d, ny = ny, 
                            x50 = x50, b4 = b4), b4 = 1,
                 ny = k - 1, k = 2)

solveE(x50, b, ny = k - 1, k = 2)

solveB(x = x50, a = 0.1, d = 2.4, x50 = 100, b4 = 1, ny = k - 1, k = 2)

}
\arguments{
    \item{x}{ a numeric vector of values at which to evaluate the function. }
    \item{y}{ a numeric vector of values at which to evaluate the function. }
    \item{a}{ a numeric parameter representing the horizontal asymptote on the
              left side (very small values of input) for \code{b} positive,
              else the horizontal asymptote on the right side:
              Zero concentration responce. }
    \item{d}{ a numeric parameter representing the horizontal asymptote on the
              right side (very large values of input) for \code{b} positive,
              else the horizontal asymptote on the left side:
              Infinite concentration responce. }
    \item{b}{ a numeric (reciprocal) scale parameter on the input axis,
              'growth rate': slope factor.
     (The 'growth' parameter \code{b} should be negative, and \code{a}
     is thus the right side (larger) asymptote, for Bertalenffy models.)
              }
  \item{x50}{ a numeric parameter representing the input value at
              the center of the curve: \eqn{IC_{50}}{IC_50}.
              The value of the \code{Richards} function will be midway
              between \code{a} and \code{d} at \code{x50}. }
    \item{e}{ a parameter determine the input value at
              the inflection point of the curve. }
   \item{ny}{ a numeric shape parameter that affects near
              which asymptote maximum 'growth' occurs.  }
    \item{k}{ an alternative numeric parameter for \code{ny}.
              (The \code{m} of Richards (1959) is \code{k}
              with \code{ny = k - 1}.) }
   \item{b4}{ a numeric parameter giving the (reciprocal) scale parameter of
              the four parameter logistic curve with slope as the 'richards'
              function at \code{x50}, and same asymptotes. }
}
\details{ 
\code{richards(x)} computes \code{d + (a - d)/z^(1/ny)} with 
\code{z = (1 + ny * (x/e)^b) = (1 + (2^ny-1) * (x/x50)^b)}
for \code{ny} different from 0.

For \code{ny} equal to 0 the limiting Gompertz curve
\code{d + (a - d) * exp(-z)} with \code{z = (x/e)^b} is found.

\code{solveE} converts the parameter \code{x50} into the parameter \code{e},
the inflection point of the curve: \code{e = x50/((2^(ny) - 1)/ny)^(1/b)}.
By \code{(2^ny-1) * (x/x50)^b = ny * (x/e)^b = ny * exp(b * (log(x)-log(e)))}
we recognize the notation of, e.g., Yin et al., 2003, and also see, how
these parameters of the above Richards function relates to the parameters
'xmid' (\code{log(x50)}) and 'scal' (\code{1/b}) of \link{SSfpl}
for \code{ny} equal 1 and the argument input the logarithm of \code{x},
\code{log(x)}.

\code{solveB} finds the parameter \code{b} for the Richard curve such
that the slope of the Richard curve is the same as the slope of
the four parameter logistic curve with parameter \code{b4} at \code{x}
(where the asymptotes \code{a} and \code{d} (and location \code{x50})
of the two curves are identical).

See \link{SSny1} for some examples on the same value found by
different functions.

We observe that to compute \code{z^(1/ny) = (1 + (2^ny-1)*(x/x50)^b)^(1/ny)}
for some values of \code{ny} the mantis \code{z = 1 + (2^ny-1)*(x/x50)^b}
should be positive.  For \code{ny > 0} this is always the case
(since \code{x} and \code{x50} are positive).
For \code{ny < 0} the power function is undefined for \code{z} negative,
unless \code{1/ny} is an integer.
The here given Richards function will for \code{ny < 0} and \code{1/ny}
an even integer have a local extreme at the asymptote \code{d},
the right hand side asymptote (\code{b} positive),
and will for \code{ny < 0} and \code{1/ny} an odd integer
just pass \code{d} for infinity.
}
% In the demo \code{sourceModified} some \code{selfStart} methods
% are presented with the functions constant the upper asymptote for
% \code{z^(1/ny) = (1 + (2^ny-1)*(x/x50)^b)^(1/ny)} negative.
\value{ A vector with function values. }
\author{ Jens Henrik Badsberg }
\references{
Richards, F.J. 1959: A flexible growth function for empirical use.
J. Exp. Bot. 10: 290-300.

Xinyou Yin, Jan Goudriaan, Egbert A. Lantinga, Jan Vos And Huub J. Spiertz:
A Flexible Sigmoid Function of Determinate Growth.
2003; Annals of Botany 91: 361-371.

Gottschalk, Paul G., Dunn, John R.:
The five-parameter logistic:
a characterization and comparison with the four-parameter logistic.
2005; Anal Biochem. Aug 1;343:54-65.
}
\keyword{ package }
\examples{
demo(showNy)

par(mfrow = c(1, 2))

# pdf("Relevant.pdf")

richardsLines
backFitFpl
relevant()

derivatives
derivatives(xlim = c(0.5, 20000), ylim = c(0, 0.1), ylim.x = c(0, 1.2), f = 4)

par(mfrow = c(2, 4))

# pdf("All.pdf")

# pdf("BelowE.pdf")
belowE(outer = FALSE, line = 1)

# pdf("Extra.pdf")
extra(outer = FALSE, line = 1)

# pdf("Below.pdf")
below(outer = FALSE, line = 1)

# pdf("Above.pdf")
above(outer = FALSE, line = 1)

for (i in 2:16) print(richards(10000, ny = -1/i))
}

