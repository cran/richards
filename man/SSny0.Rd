\name{SSny0}

% \alias{SSny.m16}
% \alias{SSny.m8}
% \alias{SSny.m4}
% \alias{SSny.m2}
  \alias{SSnym1}      % Monomolecular curve, Bertalenffy (length)
% \alias{SSny.m4o5}
% \alias{SSny.m3o4}
% \alias{SSny.m2o3}
% \alias{SSny.m1o2}
  \alias{SSnym1o3}    % Bertalenffy (weight)
% \alias{SSny.m1o4}
% \alias{SSny.m1o5}
% \alias{SSny.m1o6}
% \alias{SSny.m1o7}
% \alias{SSny.m1o8}
  \alias{SSny0}       % Gompertz
  \alias{SSny0Log}
  \alias{SSny0LogB}
% \alias{SSny.1o8}
% \alias{SSny.1o7}
% \alias{SSny.1o6}
% \alias{SSny.1o5}
% \alias{SSny.1o4}
% \alias{SSny.1o3}
% \alias{SSny.1o2}
  \alias{SSny1}       % Autocatalytic or logistic (fpl)
  \alias{SSny2}
  \alias{SSny3}
  \alias{SSny4}
% \alias{SSny.16}
% \alias{SSny.32}
% \alias{SSny.64}
% \alias{SSny.128}

\title{ 
Richards Function Models with specific values of the shape parameter 'ny' }
\description{
This \code{selfStart} model evaluates a Richards function and its gradient. 
It has an initial attribute that will evaluate initial estimates of the
parameters \code{a}, \code{d}, \code{b}, and \code{x50} for a given set
of data for \code{ny} fixed to -1, -1/3, 0, 1, 2, 3, or 4.

For \code{ny} equal to \code{-1} we have the Monomolecular curve,
also known as the the Bertalenffy model (for the length).
For \code{ny} equal to \code{-1/3} we the Bertalenffy model (for the weight).
(The 'growth' parameter \code{b} should be negative, and \code{a}
is thus the right side (larger) asymptote for Bertalenffy models.)
The Gompertz growth model is achieved by \code{ny} equal to \code{0},
and the four parameter logistic curve by  \code{ny} equal to \code{1}.
}

\usage{
   SSnym1(input, a, d,    b,  x50)

 SSnym1o3(input, a, d,    b,  x50)

    SSny0(input, a, d,    b,  x50)
 SSny0Log(input, a, d, xmid, scal)
SSny0LogB(input, a, d, xmid,    b)

    SSny1(input, a, d,    b,  x50)
    SSny2(input, a, d,    b,  x50)
    SSny3(input, a, d,    b,  x50)
    SSny4(input, a, d,    b,  x50)
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
               \code{scal} for the four point logistic curve by \code{SSfpl}. }
   \item{x50}{ a numeric parameter representing the input value at
               the center of the curve:
               The value of the function will be midway
               between \code{a} and \code{d} at \code{x50}. }
  \item{xmid}{ The value of \code{SSny0Log} will be midway between
               \code{a} and \code{d} at \code{xmid}. }
  \item{scal}{ a numeric scale parameter on the input axis, the 'growth rate'. }
}
\value{
a numeric vector of the same length as \code{input}. 
It is the value of the expression
\tabular{ll}{
 \code{ny = -1}:
    \tab \code{d + (a - d) * (1 - 1/2 * (input/x50)^b)}, \cr
    \tab if \code{1 - 1/2 * (input/x50)^b > 0},  \cr
    \tab else \code{d}, \cr
 \code{ny = -1/3}:
    \tab \code{d + (a - d) * (1 + (2^(-1/3)-1) * (input/x50)^b)^3}, \cr
    \tab if \code{1 + (2^(-1/3)-1) * (input/x50)^b > 0}, \cr
    \tab else \code{d}, \cr
 \code{ny = 0}, Log:
    \tab \code{d + (a - d) * exp( -log(2) *  exp((input - xmid) / scal))}, \cr
 \code{ny = 0}, LogB:
    \tab \code{d + (a - d) * exp(- log(2) *  exp(b * (input - xmid)))}, \cr
 \code{ny = 0}:
    \tab \code{d + (a - d) * exp( -log(2) * (input/x50)^b)}, \cr
 \code{ny = 1}:
    \tab \code{d + (a - d) / (1 + (input/x50)^b)}, \cr
 \code{ny = 2}: \tab \code{d + (a - d) / sqrt(1 + 3 * (input/x50)^b)}, \cr
 \code{ny = 3}: \tab \code{d + (a - d) / (1 + 7 * (input/x50)^b)^(1/3)}, \cr
 \code{ny = 4}: \tab \code{d + (a - d) / sqrt(sqrt((1 + 15 * (input/x50)^b)))}.
}
If all of the arguments \code{a}, \code{d}, \code{b} (\code{scal}),
and \code{x50} (\code{xmid}) are names of objects,
the gradient matrix with respect to these
names is attached as an attribute named gradient.
}
\author{ Jens Henrik Badsberg }
\keyword{ package }
\examples{ 
   x <- c(25, 50, 100, 200, 400)

            fpl(        x, 0.1, 2.4,  100,      2)
          SSny1(        x, 0.1, 2.4,  2,      100)
      SSnyFixed(1,      x, 0.1, 2.4,  2,      100)
          SSfpl(   log(x), 0.1, 2.4,  log(100), 1/2)
       richards(        x, 0.1, 2.4,  x50 = 100, b = 2, ny = 1)
          SSfpl(  -log(x), 2.4, 0.1, -log(100), 1/2)
  SSrichardsLog(   log(x), 0.1, 2.4,  log(100), 1/2, 1)
   SSnyFixedLog(1, log(x), 0.1, 2.4,  log(100), 1/2)

   a <- 2.5
   d <- 0
 x50 <- 100
   b <- 2
   e <- solveE(x50 = x50, b = b, ny = 0)
Asym <- 2.5
  b2 <- 1/e^b
  b3 <- exp(b)

   richards(     x, 2.5, 0.0, x50 = 100, b = 2, ny = 0)
      SSny0(     x, 2.5, 0.0, 2, 100)
 SSgompertz(log(x), Asym, b2, b3)
   SSny0Log(log(x), 2.5,  0.0, log(100), 1/2)


# The following calls does return identical values:

   x <- c(25, 50, 100, 200, 400)

       richards(        x, 0.1, 2.4, x50 = 100, b = 2, ny = 2)
  SSrichardsLog(   log(x), 0.1, 2.4,   log(100),  1/2, ny = 2)


# The following calls does return identical values:

   x <- c(25, 50, 100, 1000, 10000)

        # x/x50 = (-x)/(-x50):
        c(   SSrichards(         x,    a = 0.1, d = 2.4,  
                        x50 =      100,     b =    2, ny = 2))
        c(   SSrichards(        -x,    a = 0.1, d = 2.4,  
                        x50 =     -100,     b =    2, ny = 2))

        # b * (log(x) - log(xmid)) = -b * (- log(x) - -log(xmid)), 
        # note the sign on 'scal':
        c(SSrichardsLog(     log(x),   a = 0.1, d = 2.4, 
                        xmid =  log(100), scal =  1/2, ny = 2))
        c(SSrichardsLog(    -log(x),   a = 0.1, d = 2.4, 
                        xmid = -log(100), scal = -1/2, ny = 2))

   (2.4-c(SSrichardsLog(    (log(x)),  a = 2.4, d = 0.1, 
                        xmid =  log(100), scal =  1/2, ny = 2)))+0.1
   (2.4-c(SSrichardsLog(   -(log(x)),  a = 2.4, d = 0.1, 
                        xmid = -log(100), scal = -1/2, ny = 2)))+0.1


# The following calls does NOT return identical values:

   x <- c(25, 50, 100, 200, 400)

  richards(        x, 0.1, 2.4,  x50 = 100, b = 1, ny = -1)
SSrichards(        x, 0.1, 2.4,  x50 = 100, b = 1, ny = -1)

sSrichards <- selfStart( ~ d + (a - d) / (1 + (2^ny-1)*(input/x50)^b)^(1/ny),
  function (mCall, data, LHS) 
    initialRichards(mCall, data, LHS), c("a", "d", "b", "x50", "ny"))
sSrichards(        x, 0.1, 2.4,  x50 = 100, b = 1, ny = -1)
}
