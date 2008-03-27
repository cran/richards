
library(richards)

# source("sourceModified.R")

demo(sourceRichards)

demo(sourceModified)

# x <- c(25, 50, 100, 200, 400)
# SSrichardsModified(        x, 0.1, 2.4,  x50 = 100, b = 2, ny = -1.1)

modifiedPlot <- function(SS.old = sSrichards, 
                         SS.new = SSrichards,
                         log = "x", xlim = range(conc),
                         a = 0.1, d = 2.4, x50 = 100, b = 1, ny = -0.25,
                         conc = 50000 / 4^((0:(9 * n))/n), n = 50) {

plot(conc, c(SSny1(conc, 0.1, 2.4,  x50 = 100, b = 1)), 
     xlab = "Concentration", ylab = "Richards functions",
     xlim = xlim, ylim = c(0, 3), log = log, type = "l")

lines(conc, c(SS.new(conc, a = a, d = d, x50 = x50, b = b, ny = -2)), 
      col = "magenta")
lines(conc, c(SS.new(conc, a = a, d = d, x50 = x50, b = b, ny = -1)), 
      col = "cyan")
lines(conc, c(SS.new(conc, a = a, d = d, x50 = x50, b = b, ny = -.501)), 
      lwd = 4, col = "yellow")
lines(conc, c(SS.new(conc, a = a, d = d, x50 = x50, b = b, ny = -.5)), 
      col = "blue")
lines(conc, c(SS.new(conc, a = a, d = d, x50 = x50, b = b, ny = -.25)), 
      col = "green")
lines(conc, c(SS.new(conc, a = a, d = d, x50 = x50, b = b, ny = -.01)), 
      col = "red")

lines(conc, c(SS.old(conc, a = a, d = d, x50 = x50, b = b, ny = -2)), 
      lty = 2, col = "magenta")
lines(conc, c(SS.old(conc, a = a, d = d, x50 = x50, b = b, ny = -1)), 
      lty = 2, col = "cyan")
lines(conc, c(SS.old(conc, a = a, d = d, x50 = x50, b = b, ny = -.501)), 
      lty = 2, lwd = 4, col = "yellow")
lines(conc, c(SS.old(conc, a = a, d = d, x50 = x50, b = b, ny = -.5)), 
      lty = 2, col = "blue")
lines(conc, c(SS.old(conc, a = a, d = d, x50 = x50, b = b, ny = -.25)), 
      lty = 2, col = "green")
lines(conc, c(SS.old(conc, a = a, d = d, x50 = x50, b = b, ny = -.01)), 
      lty = 2, col = "red")

}


modifiedLogPlot <- function(SS.old = sSrichardsLog, 
                            SS.new = SSrichardsLog,
                            log = "", xlim = range(x), x = log(conc),
                            a = 0.1, d = 2.4, x50 = 100, b = 1, ny = -0.25,
                            conc = 50000 / 4^((0:(9 * n))/n), n = 50) {

plot(x, c(SSnyFixedLog(x, a = a, d = d,  xmid = log(100), scal = 1, ny = 1)), 
     xlab = "Log(Concentration)", ylab = "Richards functions",
     xlim = xlim, ylim = c(0, 3), log = log, type = "l", lwd = 2)

lines(x, c(SS.new(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -2)), 
      col = "magenta")
lines(x, c(SS.new(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -1)), 
      col = "cyan")
lines(x, c(SS.new(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.501)), 
      lwd = 4, col = "yellow")
lines(x, c(SS.new(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.5)), 
      col = "blue")
lines(x, c(SS.new(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.25)), 
      col = "green")
lines(x, c(SS.new(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.01)), 
      col = "red")

lines(x, c(SS.old(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -2)), 
      lty = 2, col = "magenta")
lines(x, c(SS.old(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -1)), 
      lty = 2, col = "cyan")
lines(x, c(SS.old(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.501)), 
      lty = 2, lwd = 4, col = "yellow")
lines(x, c(SS.old(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.5)), 
      lty = 2, col = "blue")
lines(x, c(SS.old(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.25)), 
      lty = 2, col = "green")
lines(x, c(SS.old(x, a = a, d = d, xmid = log(x50), scal = 1/b, ny = -.01)), 
      lty = 2, col = "red")

}


bothModifiedPlots <- function() {

modifiedPlot()
modifiedPlot(SS.old = sS.nyFixed, SS.new = SSnyFixed)
modifiedLogPlot(SS.old = sS.richardsLog, SS.new = SSrichardsLog)
modifiedLogPlot(SS.old = sS.nyFixedLog, SS.new = SSnyFixedLog)

}

show <- function(ny, a.1, a.2, a.3, label, rlim = 10^-13, dlim = 10^-14) {

  # print(label)

  # print(a.1)
  # print(a.2)
  # print(a.3)

  # print(c(a.1 - a.3))
  # print(c(a.2 - a.3))

  if (ny > 0) {
    d1 <- attributes(a.1)$gradient - attributes(a.3)$gradient
    s1 <- attributes(a.1)$gradient + attributes(a.3)$gradient
    r1 <- ifelse(d1 == 0, 0, d1 / s1)
    r1 <- ifelse(abs(r1) < rlim, 0, r1)
    if ((any(range(r1, na.rm = TRUE)) != 0) & (any(abs(range(d1, na.rm = TRUE)) > dlim))) {
      print(paste("1/3 -", label, ": ", paste(range(r1, na.rm = TRUE), collapse = ";"),
                                 " + ", paste(range(d1, na.rm = TRUE), collapse = ";")
           ))
    # print(abs(range(d1)))
    # print(r1)
    # print(a.1)
    # print(a.3)
    }
  }

  d2 <- attributes(a.2)$gradient - attributes(a.3)$gradient
  s2 <- attributes(a.2)$gradient + attributes(a.3)$gradient
  r2 <- ifelse(d2 == 0, 0, d2 / s2)
  r2 <- ifelse(abs(r2) < rlim, 0, r2)
  if ((any(range(r2, na.rm = TRUE)) != 0) & (any(abs(range(d2, na.rm = TRUE)) > dlim))) {
    print(paste("2/3 -", label, ": ", paste(range(r2, na.rm = TRUE), collapse = ";"),
                               " + ", paste(range(d2, na.rm = TRUE), collapse = ";")
         ))
  # print(abs(range(d2)))
  # print(r2)
  # print(a.2)
  # print(a.3)
  }
  # print(r2)

}

fewTests <-  function(a = 0.1, d = 2.4, x50 = 100, b = 2, ny = -0.25,
                      x = c(25, 50, 100, 200, 400, 800)) {

  subtitle <- paste("a =", a, ", d =",   d, ", x50 =", x50, 
                  ", b =", b, ", and ny =",  ny, ".")

  print(subtitle)

  if (ny > 0)
    s.2.1 <-      sS.richards(     x, a, d, x50 = x50, b = b, ny = ny)
  s.2.2 <-         SSrichards(     x, a, d, x50 = x50, b = b, ny = ny)
  s.2.3 <-    sSrichards.Cond(     x, a, d, x50 = x50, b = b, ny = ny)
  show(ny, s.2.1, s.2.2, s.2.3, "Richards")

  if (ny > 0)
    s.1.1 <-     sS.richardsE(     x, a, d, e = solveE(x50 = x50, b = b, ny = ny), b = b, ny = ny)
  s.1.2 <-        SSrichardsE(     x, a, d, e = solveE(x50 = x50, b = b, ny = ny), b = b, ny = ny)
  s.1.3 <-   sSrichardsE.Cond(     x, a, d, e = solveE(x50 = x50, b = b, ny = ny), b = b, ny = ny)
  show(ny, s.1.1, s.1.2, s.1.3, "RichardsE")

  if (ny > 0)
    s.4.1 <-   sS.richardsLog(log(x), a, d, xmid = log(x50), scal = 1 / b, ny = ny)
  s.4.2 <-      SSrichardsLog(log(x), a, d, xmid = log(x50), scal = 1 / b, ny = ny)
  s.4.3 <- sSrichardsLog.Cond(log(x), a, d, xmid = log(x50), scal = 1 / b, ny = ny)
  show(ny, s.4.1, s.4.2, s.4.3, "RichardsLog")

  if (ny > 0) {
  sSrichardsLogBM.Cond <- SSrichardsLogBM
  if (ny > 0)
    s.6.1 <-   sS.richardsLogBM(log(x), a, d, xmid = log(x50), b = b, m = 1 / ny)
  s.6.2 <-      SSrichardsLogBM(log(x), a, d, xmid = log(x50), b = b, m = 1 / ny)
  s.6.3 <- sSrichardsLogBM.Cond(log(x), a, d, xmid = log(x50), b = b, m = 1 / ny)
  show(ny, s.6.1, s.6.2, s.6.3, "RichardsBMLog")
  }

  if (ny > 0)
    s.3.1 <-       sS.nyFixed(     x, a, d, x50 = x50, b = b, ny = ny)
  s.3.2 <-          SSnyFixed(     x, a, d, x50 = x50, b = b, ny = ny)
  s.3.3 <-     sSnyFixed.Cond(     x, a, d, x50 = x50, b = b, ny = ny)
  show(ny, s.3.1, s.3.2, s.3.3, "nyFixed")

  if (ny > 0)
    s.5.1 <-    sS.nyFixedLog(log(x), a, d, xmid = log(x50), scal = 1 / b, ny = ny)
  s.5.2 <-       SSnyFixedLog(log(x), a, d, xmid = log(x50), scal = 1 / b, ny = ny)
  s.5.3 <-  sSnyFixedLog.Cond(log(x), a, d, xmid = log(x50), scal = 1 / b, ny = ny)
  show(ny, s.5.1, s.5.2, s.5.3, "nyFixedLog")

}

  fewTests(ny =    1)
  fewTests(ny =    0.5)
  fewTests(ny =    0.25)
  fewTests(ny =    0.05)
  fewTests(ny =    0.01)
  fewTests(ny =    0)
  fewTests()
  fewTests(ny =   -0.01)
  fewTests(ny =   -0.05)
  fewTests(ny =   -0.25)
  fewTests(ny =   -0.5)
  fewTests(ny =   -1)
  fewTests(ny =   -2)
  fewTests(ny =   -4)
  fewTests(ny =   -8)
  fewTests(ny =   -9)
  fewTests(ny =  -11)
  fewTests(ny =  -12)
  fewTests(ny =  -16)
  fewTests(ny =  -32)

  fewTests(ny =  -10)
  fewTests(ny =  -20)
  fewTests(ny =  -40)
  fewTests(ny = -100)

q()
