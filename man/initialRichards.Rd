\name{initialRichards}

\alias{initialNyFixed}
\alias{initialNyFixedLog}
\alias{initialNyFixedLogB}
\alias{initialRichards}
\alias{initialRichardsLog}
\alias{initialRichardsLogB}

\title{ Initial values for Richards curves}
\description{
  These functions returns initial values for the parameters
  in the Richards functions models.
}
\usage{

     initialNyFixed(mCall, data, LHS, ny = NULL, minnrow = 5,
                    pnames = c("a", "d", "b", "x50"))

  initialNyFixedLog(mCall, data, LHS, ny = NULL, minnrow = 5,
                    pnames = c("a", "d", "xmid", "scal"))

 initialNyFixedLogB(mCall, data, LHS, ny = NULL, minnrow = 5,
                    pnames = c("a", "d", "xmid", "b"))

    initialRichards(mCall, data, LHS, ny = 1, minnrow = 6, factor = 1.35, 
                    pnames = c("a", "d", "b", "x50", "ny"))

 initialRichardsLog(mCall, data, LHS, ny = 1, minnrow = 6, 
                    pnames = c("a", "d", "xmid", "scal", "ny"))

initialRichardsLogB(mCall, data, LHS, ny = 1, minnrow = 6, 
                    pnames = c("a", "d", "xmid", "b", "ny"))

}

\arguments{
    \item{mCall}{ See \link{selfStart}: "a matched call to the function model". }
     \item{data}{ See \link{selfStart}: "a data frame in which to interpret
                  the variables in mCall". }
      \item{LHS}{ See \link{selfStart}: "the expression from the left-hand side
                  of the model formula in the call to nls". }
       \item{ny}{ a numeric value, or \code{NULL}. }
  \item{minnrow}{ an integer, the number of observations required. }
   \item{factor}{ a numeric value, a factor for the independent variable. }
   \item{pnames}{ See \link{selfStart}:  "a character vector specifying the
                  terms on the right hand side of model for which initial
                  estimates should be calculated.
                  Passed as the namevec argument to the deriv function." }
}
\details{
 A four parameter logistic model is fitted to the data by
 \link{SSfpl} and \link{nls}.

 These fitted four parameters will together with \code{ny} set to 
 \code{1} determine the initial values found by \code{initialRichards},
 \code{initialRichardsLog} and \code{initialRichardsLogB}.

 For fixed \code{ny} the asymptotes and location parameter
 (the input value where the function value will be midway
 between the asymptotes)
 of a fitted four parameter model will be returned.
 A scale parameter
 (for the Richards curve with parameter \code{ny})
 will be determine such that
 the slope this curve
 is as
 the slope of the fitted four parameter curve,
 the slops for both curves midway between the asymptotes.
}
\value{
  Vector of initial values: 
  "This function should return initial values for the parameters in the model."
}
\author{ Jens Henrik Badsberg }
\keyword{ package }
\examples{

demo(sourceNy)

demo(sourceNyExtras)

demo(sourceRichards)

ls()
}