\name{simulatedPlates}
\alias{simulatedPlates}
\docType{data}
\title{ Simulated ELISA data }
\description{
  Simulated data (for an ELISA, Enzyme-Linked ImmunoSorbent Assay), 
  where the curve is a Richards curve wery close to a Gompertz curve.
}
\usage{data(simulatedPlates)}
\format{
  A data frame with 3584 observations on the following 10 variables.
  \describe{
    \item{\code{Id}}{a factor with levels \code{0.0.A} \code{0.0.B} \code{0.0.C} \code{0.0.D} \code{0.0.E} \code{0.0.F} \code{0.0.Std} \code{0.1.A} \code{0.1.B} \code{0.1.C} \code{0.1.D} \code{0.1.E} \code{0.1.F} \code{0.1.Std} \code{0.2.A} \code{0.2.B} \code{0.2.C} \code{0.2.D} \code{0.2.E} \code{0.2.F} \code{0.2.Std} \code{0.3.A} \code{0.3.B} \code{0.3.C} \code{0.3.D} \code{0.3.E} \code{0.3.F} \code{0.3.Std} \code{1.0.A} \code{1.0.B} \code{1.0.C} \code{1.0.D} \code{1.0.E} \code{1.0.F} \code{1.0.Std} \code{1.1.A} \code{1.1.B} \code{1.1.C} \code{1.1.D} \code{1.1.E} \code{1.1.F} \code{1.1.Std} \code{1.2.A} \code{1.2.B} \code{1.2.C} \code{1.2.D} \code{1.2.E} \code{1.2.F} \code{1.2.Std} \code{1.3.A} \code{1.3.B} \code{1.3.C} \code{1.3.D} \code{1.3.E} \code{1.3.F} \code{1.3.Std} \code{2.0.A} \code{2.0.B} \code{2.0.C} \code{2.0.D} \code{2.0.E} \code{2.0.F} \code{2.0.Std} \code{2.1.A} \code{2.1.B} \code{2.1.C} \code{2.1.D} \code{2.1.E} \code{2.1.F} \code{2.1.Std} \code{2.2.A} \code{2.2.B} \code{2.2.C} \code{2.2.D} \code{2.2.E} \code{2.2.F} \code{2.2.Std} \code{2.3.A} \code{2.3.B} \code{2.3.C} \code{2.3.D} \code{2.3.E} \code{2.3.F} \code{2.3.Std} \code{3.0.A} \code{3.0.B} \code{3.0.C} \code{3.0.D} \code{3.0.E} \code{3.0.F} \code{3.0.Std} \code{3.1.A} \code{3.1.B} \code{3.1.C} \code{3.1.D} \code{3.1.E} \code{3.1.F} \code{3.1.Std} \code{3.2.A} \code{3.2.B} \code{3.2.C} \code{3.2.D} \code{3.2.E} \code{3.2.F} \code{3.2.Std} \code{3.3.A} \code{3.3.B} \code{3.3.C} \code{3.3.D} \code{3.3.E} \code{3.3.F} \code{3.3.Std} \code{4.0.A} \code{4.0.B} \code{4.0.C} \code{4.0.D} \code{4.0.E} \code{4.0.F} \code{4.0.Std} \code{4.1.A} \code{4.1.B} \code{4.1.C} \code{4.1.D} \code{4.1.E} \code{4.1.F} \code{4.1.Std} \code{4.2.A} \code{4.2.B} \code{4.2.C} \code{4.2.D} \code{4.2.E} \code{4.2.F} \code{4.2.Std} \code{4.3.A} \code{4.3.B} \code{4.3.C} \code{4.3.D} \code{4.3.E} \code{4.3.F} \code{4.3.Std} \code{5.0.A} \code{5.0.B} \code{5.0.C} \code{5.0.D} \code{5.0.E} \code{5.0.F} \code{5.0.Std} \code{5.1.A} \code{5.1.B} \code{5.1.C} \code{5.1.D} \code{5.1.E} \code{5.1.F} \code{5.1.Std} \code{5.2.A} \code{5.2.B} \code{5.2.C} \code{5.2.D} \code{5.2.E} \code{5.2.F} \code{5.2.Std} \code{5.3.A} \code{5.3.B} \code{5.3.C} \code{5.3.D} \code{5.3.E} \code{5.3.F} \code{5.3.Std} \code{6.0.A} \code{6.0.B} \code{6.0.C} \code{6.0.D} \code{6.0.E} \code{6.0.F} \code{6.0.Std} \code{6.1.A} \code{6.1.B} \code{6.1.C} \code{6.1.D} \code{6.1.E} \code{6.1.F} \code{6.1.Std} \code{6.2.A} \code{6.2.B} \code{6.2.C} \code{6.2.D} \code{6.2.E} \code{6.2.F} \code{6.2.Std} \code{6.3.A} \code{6.3.B} \code{6.3.C} \code{6.3.D} \code{6.3.E} \code{6.3.F} \code{6.3.Std} \code{7.0.A} \code{7.0.B} \code{7.0.C} \code{7.0.D} \code{7.0.E} \code{7.0.F} \code{7.0.Std} \code{7.1.A} \code{7.1.B} \code{7.1.C} \code{7.1.D} \code{7.1.E} \code{7.1.F} \code{7.1.Std} \code{7.2.A} \code{7.2.B} \code{7.2.C} \code{7.2.D} \code{7.2.E} \code{7.2.F} \code{7.2.Std} \code{7.3.A} \code{7.3.B} \code{7.3.C} \code{7.3.D} \code{7.3.E} \code{7.3.F} \code{7.3.Std}}
    \item{\code{PT}}{Sampletype: a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{Std}}
    \item{\code{LAY}}{a numeric vector: 0, 1, 2, ..., 7}
    \item{\code{PLA}}{Plate: a numeric vector: 0, 1, 2, ..., 3}
    \item{\code{PR}}{platerow: a numeric vector: 0 or 1}
    \item{\code{DS}}{Delutionstep: a numeric vector: 0, 1, 2, ..., 7}
    \item{\code{CONC}}{Concentration: a numeric vector}
    \item{\code{LOGFORT}}{Logarithm of concentration: a numeric vector}
    \item{\code{MEAN}}{The simulated mean on the curve of the sample: a numeric vector}
    \item{\code{SIGNAL}}{The simulated absorbance, signal: a numeric vector}
  }
}
% \details{ }
% \source{ }
% \references{ }
\examples{data(simulatedPlates)}
\keyword{datasets}
