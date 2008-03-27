
# library(richards)

# source("sourceModified.R")

demo(sourceRichards)

demo(sourceModified)

toMatrix <- function(A) {
  dims <- dim(A)
  wch <- dims != 1
  M <- matrix(A, dims[wch])
  dimnames(M) <- dimnames(A)[wch]
  return(M)
}

gradArray <- function(var = "a", logx = TRUE, xlim = range(conc),
                      SS = SSrichards,
                      a = 0.1, d = 2.4, x50 = 100, b = 2, ny = 1,
                      conc = 100 * 10^seq(-2, 2, length = 1 + 2 * n), n = 25,
                      l = length(conc)) {
  dims <- c(l, length(a), length(d), length(x50), length(b), length(ny))
  dimn <- list(paste(format(ifelse(rep(logx, length(conc)), log10(conc), conc),
                            scientific = FALSE, digits = 5)),
               paste(a), paste(d), paste(x50), paste(b), paste(ny))
  names(dims) <- c("Conc", "a", "d", "x50", "b", "ny")
  wch <- dims != 1
  wf <- names(dims)[wch]
  if (!(length(wf) == 2)) warning("To match only one variable should have more values")
  f <- function(x) rep(x, rep(l, length(x)))
  S <- SS(conc, a = f(a), d = f(d), x50 = f(x50), b = f(b), ny = f(ny))
  grad <- attributes(S)$gradient[, var]
  A <- array(grad, dims)
  dimnames(A) <- dimn
  return(A)
}

myContour <- function(var = "a", logx = TRUE, logy = TRUE, SS = SSrichards,
                      a = 0.1, d = 2.4, x50 = 100, b = 2, ny = 1,
                      conc = 100 * 10^seq(-2, 2, length = 1 + 2 * n), n = 25,
                      xaxis = FALSE) {
  A <- gradArray(var = var, logx = logx,
                 SS = SS, a = a, d = d, x50 = x50, b = b, ny = ny, conc = conc)
  dims <- dim(A)
  wch <- dims != 1
  fnames <- c("Conc", "a", "d", "x50", "b", "ny")
  wf <- fnames[wch]
  if (!(length(wf) == 2))
    stop("One variable should have more values")
  ylab <- wf[2]
  M <- toMatrix(A)
  x <- as.real(dimnames(M)[[1]])
  y <- as.real(dimnames(M)[[2]])
  if (logy & (min(y) > 0)) {
    y <- log10(y)
    ylab <- paste("Log(", ylab, ")")
  }
  if ((ylab == "ny") & logy & (max(y) < 0)) {
    y <- -log10(-y)
    ylab <- paste("-Log(-", ylab, ")")
  }
  main <- paste("Grad.", var)
  mar <- par()$mar
  if (var == "x50")
    M <- -M
  # print(c(var, range(M)))
  if ((var == "a") | (var == "d") | (var == "x50")) {
    M <- log10(M)
    main <- paste("Log(", main, ")")
  }
  if (mar[3] < 1)
    main <- ""
  if (mar[2] < 1)
    ylab <- ""
  contour(x = x, y = y, z = M, main = main, axes = FALSE,
          xlab = "Log(Concentration)", ylab = ylab)
  if (mar[2] > 1)
    axis(side = 2)
  if (xaxis)
    axis(side = 1)
  mar[2] <- 0.1
  par(mar = mar)
}

# b <-   2^-3 * sqrt(2)^((-m:m)/2)

# SS   <-  SSrichards
# A    <-  0.1
# D    <-  2.4
# X50  <-  100
# B    <-  2
# Ny   <-  1
# m    <-  24
# n    <-  50/9
# conc <-  100 * 10^seq(-2, 2, length = 1 + 2 * n)

contour5by5 <- function(SS = SSrichards, title = "SS?",
                     A = 0.1, D = 2.4, X50 = 100, B = 2, Ny = 1, m = 48,
                     conc = 100 * 10^seq(-2, 2, length = 1 + 2 * n), n = 25) {

  subtitle <- paste("a =", A, ", d =",   D, ", x50 =", X50, 
                  ", b =", B, ", and ny =",  Ny, ".")
  # subtitle <<- subtitle
  a    <-     A + seq(-0.25, 0.25, length = 2 * m)
  d    <-     D + seq(-0.25, 0.25, length = 2 * m)
  x50  <-  rev(X50 / 100 * 50000 / 4^(9 * seq(0, 1, length = 2 * m)))
  x50  <-  rev(X50 / 100 * 50000 / 4^(9 * seq(0, 1, length = 2 * m)))
  b    <-     B * 10^seq(-2, 2, length = 2 * m)
  ny   <-    Ny *  2^(-m:m)
  ny   <-    Ny * 10^seq(-4, 2, length = 2 * m)
  a    <- sort(a)
  d    <- sort(d)
  x50  <- sort(x50)
  b    <- sort(b)
  ny   <- sort(ny)
  vars <- c("a", "x50", "b", "ny", "d")
  par(mfrow = c(5, 5))
  par(oma = c(1, 0, 2, 0))
  par(mar = c(0.1, 4.1, 2.1, 0.1))
  for (var in vars) 
    myContour(SS = SS, var = var, a = a, logy = FALSE)
  par(mar = c(0.1, 4.1, 0.1, 0.1))
  for (var in vars) 
    myContour(SS = SS, var = var, x50 = x50, logy = TRUE)
  par(mar = c(0.1, 4.1, 0.1, 0.1))
  for (var in vars) 
    myContour(SS = SS, var = var, b = b, logy = TRUE)
  par(mar = c(0.1, 4.1, 0.1, 0.1))
  for (var in vars) 
    myContour(SS = SS, var = var, ny = ny, logy = TRUE)
  par(mar = c(2.1, 4.1, 0.1, 0.1))
  for (var in vars) 
    myContour(SS = SS, var = var, d = d, logy = FALSE, xaxis = TRUE)
  title(main = title, outer = TRUE)
  title(sub = subtitle, outer = TRUE, line = 0)
}
  
make3plots <- function() {
  pdf("contours.pdf")
  contour5by5(title = "SSrichards, Ny = 1")
  contour5by5(Ny = -0.1, title = "SSrichards, Ny = -0.1")
  contour5by5(SS = sSrichards, Ny = -0.1, title = "sSrichards, Ny = -0.1")
}

# make3plots()

