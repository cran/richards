\name{SSrichards}
\alias{SSrichards}
\alias{SSrichardsE}
\alias{SSrichardsLog}
\alias{SSrichardsLogBM}
\alias{SSrichards}
\title{ Richards Function Model }
\description{
This selfStart model evaluates the Richards function and its gradient,
for \code{ny} different from \code{0}. 
It has an initial attribute that will evaluate initial estimates
of the parameters \code{a}, \code{d}, \code{b}, and \code{x50}
for a given set of data for \code{ny} set to \code{1}. 
}
\usage{
             SSrichards(input, a, d, b,     x50, ny)
            SSrichardsE(input, a, d, b,       e, ny)
          SSrichardsLog(input, a, d, xmid, scal, ny)
        SSrichardsLogBM(input, a, d, xmid,    b,  m)
}
\arguments{
  \item{input}{ a numeric vector of values at which to evaluate the model. }
      \item{a}{ a numeric parameter representing the horizontal asymptote on the
                left side (very small values of input) for \code{b} positive,
                else the horizontal asymptote on the right side. }
      \item{d}{ a numeric parameter representing the horizontal asymptote on the
                right side (very large values of input) for \code{b} positive,
                else the horizontal asymptote on the left side. }
      \item{b}{ a numeric scale parameter on the input axis,
                the 'growth rate', the reciprocal of the scale parameter
                \code{scal} for the four point logistic curve by SSfpl. }
    \item{x50}{ a numeric parameter representing the input value at
                the inflection point (for \code{ny = 1}) of the curve.
                The value of \code{SSrichards} will be
                midway between \code{a} and \code{d} at \code{x50}. }
      \item{e}{ the parameter determine the input value at
                the inflection point of the curve. }
     \item{ny}{ the fifth parameter: affects near which asymptote maximum
                'growth' occurs. }
      \item{m}{ the fifth parameter: affects near which asymptote maximum
                'growth' occurs. Reciprocal of \code{ny}.
                (Note that the \code{m} of Richards (1959) is \code{k}
                with \code{ny = k - 1}.) }
   \item{xmid}{ a numeric parameter representing the input value at the 
                inflection point (for \code{ny = 1}) of the curve.
                The value of \code{SSrichardsLog} will be
                midway between \code{a} and \code{d} at \code{xmid}. }
   \item{scal}{ a numeric scale parameter on the input axis, 'growth rate'. }
}
\value{
a numeric vector of the same length as \code{input}. 
It is the value of the expression
\tabular{ll}{
% SSrichards:      \tab \code{d + (a - d) / (1 + (2^ny-1) * (input/x50)^b)^(1/ny)}, \cr
% SSrichardsE:     \tab \code{d + (a - d) / (1 +       ny * (input/e)^b)^(1/ny)}, \cr
% SSrichardsLog:   \tab \code{d + (a - d) / (1 + (2^ny-1) * exp((input-xmid)/scal))^(1/ny)}, \cr
% SSrichardsLogBM: \tab \code{d + (a - d) / (1 + (2^(1/m)-1) * exp(b*(input-xmid)))^m}, \cr
SSrichards:        \tab \code{d + (a - d) / (1 + (2^ny-1) * (input/x50)^b)^(1/ny)} \cr
                   \tab if \code{1 + (2^ny-1) * (input/x50)^b > 0},  \cr
                   \tab else \code{d}, \cr
SSrichardsE:       \tab \code{d + (a - d) / (1 +       ny * (input/e)^b)^(1/ny)}, \cr
                   \tab if \code{1 +       ny * (input/e)^b > 0}, \cr
                   \tab else \code{d}, \cr
SSrichardsLog:     \tab \code{d + (a - d) / (1 + (2^ny-1) * exp((input-xmid)/scal))^(1/ny)} \cr
                   \tab if \code{1 + (2^ny-1) * exp((input-xmid)/scal) > 0}, \cr
                   \tab else \code{d}, or \cr
SSrichardsLogBM:   \tab \code{d + (a - d) / (1 + (2^(1/m)-1) * exp(b*(input-xmid)))^m} \cr
                   \tab if \code{1 + (2^(1/m)-1) * exp(b*(input-xmid)) > 0}, \cr
                   \tab else \code{d}. 
}

See \link{richards} about details of the function.
Observe that the \code{selfStarts} models \code{SSrichards}, \code{SSrichardsE},
and \code{SSrichardsLog} does not handle the Gompertz function.
See also \link{richards} about the modified \code{selfStarts} with constant 
asymptote on the right side (\code{b} positive) after reaching this asymptote.

If all of the arguments \code{a}, \code{d}, \code{b} (\code{scal}),
\code{x50} (\code{xmid}), and \code{ny} are names of objects, the gradient
matrix with respect to these names is attached as an attribute named gradient.
}
\details{
Why the 4 versions?

The version \code{SSrichardsE} (and mainly \code{SSrichards}) was made to be
close to the usual parameterization of the four-parameter logistic curve
at Statens Serum Institute. 

\code{SSrichardsLog} and \code{SSrichardsLogBM} was made during writing of this
documentation.
\code{SSrichardsLog} is closer to other \code{selfStart} methods implemented in R.
For \code{SSrichardsLogBM} bounds on  the parameter \code{b} can be avioded.
}
\author{ Jens Henrik Badsberg }
\keyword{ package }
\examples{ 

demo(contourModified)



# make3plots()
}
