
# library(richards)

`backFitFpl` <-
function (k, conc = z, z = 100 * 2^(-7:7), y = richards(z, k = k),
          a = 2.4, d = 0.1, x50 = 100, b4 = -1) 
{
    # OK:
    fit <- nls(y ~ SSny1(x, a, d, b, x50),
               start = list(a = a, d = d, b = b4, x50 = x50),
               data = data.frame(y = y, x = z),
               nls.control(maxiter = 100))
    coef <- summary(fit)$coefficients[, "Estimate"]
    backfitted <- fpl.inv(richards(conc, k = k), a = coef["a"], 
                   d = coef["d"], e = coef["x50"], b = coef["b"])
    return(backfitted)
}

`richardsLines` <-
function (i, kv, conc, col, k = kv[i], l = length(kv), cex = 0.75,
          a = 2.4, d = 0.1, x50 = 100, b4 = -1) 
{
    setLabel <- function(k, x, txt, j) {
        p <- -round(log10(abs(x))) + 2
        if (x != 0) 
            x <- round(x * 10^p)/10^p
        text(ifelse(k > 130, 3^j, max(conc)/3^j), 2.5 * (l - 0)/(l * 2), 
             txt, col = "black", cex = cex)
        text(ifelse(k > 130, 3^j, max(conc)/3^j), 2.5 * (l - i)/(l * 2), 
             paste(x), col = col, cex = cex)
    }
    # print(k)
    setLabel(k, k,         "K",  1)
    setLabel(k, k - 1,     "Ny", 4)
  # setLabel(k, 2 - k,     "l",  0)
  # setLabel(k, 1/(k - 1), "m",  3)

    ny <- k - 1
    b <- solveB(ny = ny, a = a, d = d, x50 = x50, b4 = b4)
    e <- solveE(x50 = 100, b = ny, ny = ny)
    y <- richards(e, k = kv[i], a = a, d = d, x50 = x50, b4 = b4)
    points(e, y, col = col, cex = 1, pch = 20)

    x0 <- 100 * (1-2^(ny))^(-1/b)
    x0 <- c(x0 / 1.25, x0, x0 * 1.25, 40000)
    y0 <- richards(x0, k = kv[i], a = a, d = d, x50 = x50, b4 = b4)
  # print(paste(format(c(k, ny, x0, y0, 1/ny), 
  #                    digits = 1, scientific = FALSE), collapse = ", "))

    points(x0[2], y0[2], col = col, cex = 1, pch = 3)

    Conc <- sort(c(conc, x0))
    lines(Conc, richards(Conc, k = k, a = a, d = d, x50 = x50, b4 = b4),
          lty = i, col = col)
}


`derivatives` <-
function (kv = NULL, ylim = c(0.0, 0.01), ylim.x = c(0.3, 0.7),
          z = 100 * 2^(-7:7), xlim = c(50, 200),
          conc = 20000 / 4^((0:(9 * n))/n), n = 50, # 401 points
          offset = 546, f = 8, cex = 0.75,
          a = 2.4, d = 0.1, x50 = 100, b4 = -1) 
{
    if (is.null(kv)) {
        kv <- c(-7, -3, -1, 0, 1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5, 5/6, 1,
                7/6, 6/5, 5/4, 4/3, 3/2, 3, 5, 9, 17, 33, 65, 129)
        kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    }

    # 1. derivarive, with respect to 'x':
    plot(conc, fpl.deriv(conc, a = a, d = d, e = x50, b = b4),
         type = "l", log = "x", xlim = xlim, xlab = "Concentration", 
         ylim = ylim, ylab = "Derivative of Richards function")
    for (i in 1:length(kv)) {
        # if (kv[i] == 0.5) 
        #   kv[i] <- 0.50000001
        col <- colors()[i * f + offset]
        lines(conc, richards.deriv(conc, k = kv[i],
                                   a = a, d = d, x50 = x50, b4 = b4),
              lty = i, col = col, cex = cex)
    }

    # 1. derivarive, with respect to 'log(x)':
    plot(conc,  fpl.deriv.x(conc, a = a, d = d, e = x50, b = b4),
         type = "l", log = "x", xlim = xlim, xlab = "Concentration", 
         ylim = ylim.x, ylab = "Derivative with resp. to log(Concentration)")
    points(100, fpl.deriv.x(100))
    for (i in 1:length(kv)) {
        # if (kv[i] == 0.5) 
        #   kv[i] <- 0.50000001
        col <- colors()[i * f + offset]
        lines(conc, richards.deriv.x(conc, k = kv[i],
                                     a = a, d = d, x50 = x50, b4 = b4),
              lty = i, col = col, cex = cex)
        b <- solveB(a = a, d = d, x50 = x50, b4 = b4, ny = kv[i] - 1)
        e <- solveE(x50 = x50, b = b, ny = kv[i] -1)
        y <- richards.deriv.x(e, a = a, d = d, b = b, e = e, k = kv[i])
        points(e, y, col = col, cex = 2, pch = 4)
    }

}

`funcBack` <-
function (kv = NULL, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7), 
          conc = 20000 / 4^((0:(9 * n))/n), n = 50, # 401 points
          offset = 546, f = 8, cex = 0.75,
          a = 2.4, d = 0.1, x50 = 100, b4 = -1) 
{
    if (is.null(kv)) {
        kv <- c(-7, -3, -1, 0, 1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5, 5/6, 1,
                7/6, 6/5, 5/4, 4/3, 3/2, 3, 5, 9, 17, 33, 65, 129)
        kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    }

    plot(conc, fpl(conc, a = a, d = d, e = x50, b = b4),
         type = "l", log = "x", lwd = 2,
         xlim = c(0.5, max(conc)), ylim = c(0, 2.5),
         xlab = "Concentration", ylab = "Richards function")

    for (i in 1:length(kv)) {
        # if (kv[i] == 0.5) 
        #    kv[i] <- 0.50000001
        col <- colors()[i * f + offset]
        richardsLines(i, kv, conc, col, cex = cex,
                      a = a, d = d, x50 = x50, b4 = b4)
    }

    plot(conc, rep(1, length(conc)), type = "l", log = "xy", 
         xlim = c(0.5, max(conc)), ylim = ylim, xlab = "Concentration", 
         ylab = "Backfitted relative concentration")

    for (i in 1:length(kv)) {
        col <- colors()[i * f + offset]
        k <- kv[i]
        p <- 2 - round(log10(abs(k)))
        k <- round(k * 10^p)/10^p
        # print(k)
        if ((kv[i] != 2) & (kv[i] > 0)) {
            bf <- rep(1, length(conc))
            # if (kv[i] == 0.5) 
            #    kv[i] <- 0.50000001
            try(bf <- backFitFpl(k = kv[i], conc, z,
                                 a = a, d = d, x50 = x50, b4 = b4)/conc)
            delta <- 2 * ((i - 1)/length(kv))^2
            idx <- (1:length(conc))[bf < 1 + delta]
            upp <- max(idx[!is.na(idx)])
            text(conc[upp], bf[upp], paste(k), col = col, cex = cex)
            low <- min(idx[!is.na(idx)])
            text(conc[low], bf[low], paste(k), col = col, cex = cex)
            lines(conc, bf, lty = i, col = col)
        }
    }
    # print(cbind(k = kv, ny = kv - 1, factor = 2^(kv - 1)-1, 
    #             potens = 1/(kv - 1), root = (kv - 1)))
}

`relevant` <- function (outer = TRUE, line = -2, cex = 0.75,
                        a = 2.4, d = 0.1, x50 = 100, b4 = -1) {
    kv <- c(-7, -3, -1, 0, 1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5, 5/6, 1,
            7/6, 6/5, 5/4, 4/3, 3/2, 3, 5, 9, 17, 33, 65, 129)
    kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    funcBack(kv = kv, cex = cex,
             ylim = c(0.6, 1.5), z = 100 * 2^(-7:7), offset = 546, f = 4,
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Ny from -16 to 128", outer = outer, line = line)
}

`extraPlus` <- function (outer = TRUE, line = -2, cex = 0.75,
                         a = 2.4, d = 0.1, x50 = 100, b4 = -1) {
    kv = c(-3, -2, -1, -7/8, -4/5, -3/4, -2/3, -1/c(2:8, 16),
                    0, 1/c(16, 8:2), 1)
    kv <- ifelse(kv == -1, 0.0, kv + 1.000000)
    funcBack(kv =  kv, f = 4,
             cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7),
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Below Fpl : Ny from -3 to 1", outer = outer, line = line)
}

`Bertalenffy` <- function (outer = TRUE, line = -2, cex = 0.75,
                           a = 2.4, d = 0.1, x50 = 100, b4 = -1) {
    kv = c( -1, -1/3)
    kv <- ifelse(kv == -1, 0.0, kv + 1.000000)
    funcBack(kv =  kv, f = 4,
             cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7),
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Bertalenffy : Ny -1 and -1/3", outer = outer, line = line)
}

`extra` <- function (outer = TRUE, line = -2, cex = 0.75,
                     a = 2.4, d = 0.1, x50 = 100, b4 = -1) {
    kv = c(-16, -8, -4, -2, -1, -4/5, -3/4, -2/3, -0.5, 0, 0.5, 1)
    kv <- ifelse(kv == -1, 0.0, kv + 1.0000001)
    funcBack(kv = kv, cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7),
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Below Fpl : Ny from -16 to 1", outer = outer, line = line)
}

`above` <- function (outer = TRUE, line = -2, cex = 0.75,
                     a = 2.4, d = 0.1, x50 = 100, b4 = -1) 
{
  # funcBack(kv = 2^((1:14)/2)+1,
  #          cex = cex, ylim = c(0.25, 2.0), z = 100 * 2^(-6:6))
    funcBack(kv = 2^((1:14)/2)+1,
             cex = cex, ylim = c(0.25, 2.0), z = 100 * 2^(-7:7),
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Above Fpl: Ny from 1 to 128", outer = outer, line = line)
}

`belowE` <- function (outer = TRUE, line = -2, cex = 0.75,
                      a = 2.4, d = 0.1, x50 = 100, b4 = -1) {
    # funcBack(kv = 2^((-8:2)/2), cex = cex, ylim = c(0.6, 1.25))
    kv <- 2^((-7:4)/4)
    kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    funcBack(kv = kv, cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7),
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Below Fpl : Ny from -0.7 to 1", outer = outer, line = line)
}

`below` <- function (outer = TRUE, line = -2, cex = 0.75,
                     a = 2.4, d = 0.1, x50 = 100, b4 = -1) {
    kv <- c(-1/(2:6)+0.000000001, 0, 1/(6:2)+0.000000001) + 1
    kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    funcBack(kv = kv,
             cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7),
             a = a, d = d, x50 = x50, b4 = b4)
    title(main = "Below Fpl : Ny from -0.5 to 0.5", outer = outer, line = line)
}

`pdfNy` <- function () {
    pdf("Relevant.pdf")
    par(mfrow = c(1, 2))
    relevant(cex = 0.5)
    derivatives(f = 4)
    derivatives(xlim = c(0.5, 20000),
                ylim = c(0, 0.1), ylim.x = c(0, 1.2), f = 4)
    relevant(cex = 0.5, a = 0.1, d = 2.4, x50 = 100, b4 = 1)
    derivatives(f = 4, a = 0.1, d = 2.4, x50 = 100, b4 = 1)
    derivatives(xlim = c(0.5, 20000),
                ylim = c(0, 0.1), ylim.x = c(0, 1.2), 
                a = 0.1, d = 2.4, x50 = 100, b4 = 1, f = 4)
    pdf("All.pdf")
    par(mfrow = c(2, 4))
    belowE(outer = FALSE, line = 1, cex = 0.5)
     extra(outer = FALSE, line = 1, cex = 0.5)
     below(outer = FALSE, line = 1, cex = 0.5)
     above(outer = FALSE, line = 1, cex = 0.5)
}


par(mfrow = c(1, 2))

# pdf("Relevant.pdf")

relevant()

derivatives(xlim = c(0.5, 20000), ylim = c(0, 0.1), ylim.x = c(0, 1.2), f = 4)

relevant(a = 0.1, d = 2.4, x50 = 100, b4 = 1)

derivatives(xlim = c(0.5, 20000), ylim = c(0, 0.1), ylim.x = c(0, 1.2), 
            a = 0.1, d = 2.4, x50 = 100, b4 = 1, f = 4)


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

# q()
