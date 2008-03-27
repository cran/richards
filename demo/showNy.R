
# library(richards)

`backFitFpl` <-
function (k, conc = z, z = 100 * 2^(-7:7), y = richards(z, k = k)) 
{
    # OK:
    fit <- nls(y ~ d + (a - d)/(1 + (x/e)^b),
               data = data.frame(y = y, x = z),
               nls.control(maxiter = 100),
               start = list(a = 0.1, d = 2.4, e = 100, b = 1))
    coef <- summary(fit)$coefficients[, "Estimate"]
    backfitted <- fpl.inv(richards(conc, k = k), a = coef["a"], 
                   d = coef["d"], e = coef["e"], b = coef["b"])
    return(backfitted)
}

`backFitFplWithout` <-
function (k, conc = z, z = 100 * 2^(-7:7), y = richards(z, k = k)) 
{
    # OK:
    fit <- nls(y ~ d + (a - d)/(1 + (x/e)^b),
               data = data.frame(y = y, x = z),
               nls.control(maxiter = 100),
               start = list(a = 0.1, d = 2.4, e = 100, b = 1))
    coef <- summary(fit)$coefficients[, "Estimate"]
    backfitted <- fpl.inv(richards(conc, k = k), a = coef["a"], 
                   d = coef["d"], e = coef["e"], b = coef["b"])
    return(backfitted)
}

`backFitFplWrong` <-
function (k, conc = z, z = 100 * 2^(-7:7), y = richards(z, k = k)) 
{
    # This does not work just because of the switch of 'b' and 'x50':
    fit <- nls(y ~ SSny1(x, a, d, b, x50),
               start = list(a = 0.1, d = 2.4, x50 = 100, b = 1),
               data = data.frame(y = y, x = z),
               nls.control(maxiter = 100))
    coef <- summary(fit)$coefficients[, "Estimate"]
    backfitted <- fpl.inv(richards(conc, k = k), a = coef["a"], 
                   d = coef["d"], e = coef["x50"], b = coef["b"])
    return(backfitted)
}

`backFitFplNS` <-
function (k, conc = z, z = 100 * 2^(-7:7), y = richards(z, k = k)) 
{
    # No starting values for 'SSny1' ?
    fit <- nls(y ~ SSny1(x, a, d, b, x50),
               data = data.frame(y = y, x = z), algorithm = "port",
               nls.control(maxiter = 50))
    coef <- summary(fit)$coefficients[, "Estimate"]
    backfitted <- fpl.inv(richards(conc, k = k), a = coef["a"], 
                   d = coef["d"], e = coef["x50"], b = coef["b"])
    return(backfitted)
}

`backFitFpl` <-
function (k, conc = z, z = 100 * 2^(-7:7), y = richards(z, k = k)) 
{
    # OK:
    fit <- nls(y ~ SSny1(x, a, d, b, x50),
               start = list(a = 0.1, d = 2.4, b = 1, x50 = 100),
               data = data.frame(y = y, x = z),
               nls.control(maxiter = 100))
    coef <- summary(fit)$coefficients[, "Estimate"]
    backfitted <- fpl.inv(richards(conc, k = k), a = coef["a"], 
                   d = coef["d"], e = coef["x50"], b = coef["b"])
    return(backfitted)
}

`richardsLines` <-
function (i, kv, conc, col, k = kv[i], l = length(kv), cex = 0.75) 
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
    lines(conc, richards(conc, k = k), lty = i, col = col)
    e <- solveE(x50 = 100, b = solveB(ny = kv[i] - 1), ny = kv[i] -1)
    y <- richards(e, k = kv[i])
    points(e, y, col = col, cex = 1, pch = 20)
}

`derivatives` <-
function (kv = NULL, ylim = c(0.0, 0.01), ylim.x = c(0.3, 0.7),
          z = 100 * 2^(-7:7), xlim = c(50, 200),
          conc = 20000 / 4^((0:(9 * n))/n), n = 50, # 401 points
          offset = 546, f = 8, cex = 0.75) 
{
    if (is.null(kv)) {
        kv <- c(-7, -3, -1, 0, 1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5, 5/6, 1,
                7/6, 6/5, 5/4, 4/3, 3/2, 3, 5, 9, 17, 33, 65, 129)
        kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    }

    # 1. derivarive, with respect to 'x':
    plot(conc, fpl.deriv(conc), type = "l", log = "x",
         xlim = xlim, xlab = "Concentration", 
         ylim = ylim, ylab = "Derivative of Richards function")
    for (i in 1:length(kv)) {
        # if (kv[i] == 0.5) 
        #   kv[i] <- 0.50000001
        col <- colors()[i * f + offset]
        lines(conc, richards.deriv(conc, k = kv[i]),
              lty = i, col = col, cex = cex)
    }

    # 1. derivarive, with respect to 'log(x)':
    plot(conc,  fpl.deriv.x(conc), type = "l", log = "x",
         xlim = xlim, xlab = "Concentration", 
         ylim = ylim.x, ylab = "Derivative with resp. to log(Concentration)")
    points(100, fpl.deriv.x(100))
    for (i in 1:length(kv)) {
        # if (kv[i] == 0.5) 
        #   kv[i] <- 0.50000001
        col <- colors()[i * f + offset]
        lines(conc, richards.deriv.x(conc, k = kv[i]),
              lty = i, col = col, cex = cex)
        b <- solveB(ny = kv[i] - 1)
        e <- solveE(x50 = 100, b = b, ny = kv[i] -1)
        y <- richards.deriv.x(e, k = kv[i])
        points(e, y, col = col, cex = 2, pch = 4)
    }

}

`funcBack` <-
function (kv = NULL, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7), 
          conc = 20000 / 4^((0:(9 * n))/n), n = 50, # 401 points
          offset = 546, f = 8, cex = 0.75) 
{
    if (is.null(kv)) {
        kv <- c(-7, -3, -1, 0, 1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5, 5/6, 1,
                7/6, 6/5, 5/4, 4/3, 3/2, 3, 5, 9, 17, 33, 65, 129)
        kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    }

    plot(conc, fpl(conc), type = "l", log = "x", lwd = 2,
         xlim = c(0.5, max(conc)), ylim = c(0, 2.5),
         xlab = "Concentration", ylab = "Richards function")

    for (i in 1:length(kv)) {
        if (kv[i] == 0.5) 
           kv[i] <- 0.50000001
        col <- colors()[i * f + offset]
        richardsLines(i, kv, conc, col, cex = cex)
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
            if (kv[i] == 0.5) 
               kv[i] <- 0.50000001
            try(bf <- backFitFpl(k = kv[i], conc, z)/conc)
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

`relevant` <- function (outer = TRUE, line = -2, cex = 0.75) {
    kv <- c(-7, -3, -1, 0, 1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5, 5/6, 1,
            7/6, 6/5, 5/4, 4/3, 3/2, 3, 5, 9, 17, 33, 65, 129)
    kv <- ifelse(kv == 1, 1, kv + 0.0000001)
    funcBack(kv = kv, cex = cex,
             ylim = c(0.6, 1.5), z = 100 * 2^(-7:7), offset = 546, f = 4)
    title(main = "Ny from -16 to 128", outer = outer, line = line)
}

`extra` <- function (outer = TRUE, line = -2, cex = 0.75) {
    funcBack(kv = c(-16, -8, -4, -2, -1, -4/5, -3/4, -2/3, -0.5, 0, 0.5, 1) + 1,
             cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7))
    title(main = "Above Fpl : Ny from -16 to 1", outer = outer, line = line)
}

`below` <- function (outer = TRUE, line = -2, cex = 0.75) 
{
  # funcBack(kv = 2^((1:14)/2)+1,
  #          cex = cex, ylim = c(0.25, 2.0), z = 100 * 2^(-6:6))
    funcBack(kv = 2^((1:14)/2)+1,
             cex = cex, ylim = c(0.25, 2.0), z = 100 * 2^(-7:7))
    title(main = "Below Fpl: Ny from 1 to 128", outer = outer, line = line)
}

`aboveE` <- function (outer = TRUE, line = -2, cex = 0.75) {
    # funcBack(kv = 2^((-8:2)/2), cex = cex, ylim = c(0.6, 1.25))
    funcBack(kv = 2^((-7:4)/4),
             cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7))
    title(main = "Above Fpl : Ny from -0.7 to 1", outer = outer, line = line)
}

`above` <- function (outer = TRUE, line = -2, cex = 0.75) {
    funcBack(kv = c(-1/(2:6)+0.000000001, 0, 1/(6:2)+0.000000001) + 1,
             cex = cex, ylim = c(0.6, 1.5), z = 100 * 2^(-7:7))
    title(main = "Above Fpl : Ny from -0.5 to 0.5", outer = outer, line = line)
}

`pdfNy` <- function () {
    pdf("Relevant.pdf")
    par(mfrow = c(1, 2))
    relevant(cex = 0.5)
    derivatives(f = 4)
    derivatives(xlim = c(0.5, 20000),
                ylim = c(0, 0.1), ylim.x = c(0, 1.2), f = 4)
    pdf("All.pdf")
    par(mfrow = c(2, 4))
    aboveE(outer = FALSE, line = 1, cex = 0.5)
     extra(outer = FALSE, line = 1, cex = 0.5)
     above(outer = FALSE, line = 1, cex = 0.5)
     below(outer = FALSE, line = 1, cex = 0.5)
}

par(mfrow = c(1, 2))

# pdf("Relevant.pdf")

relevant()


derivatives(xlim = c(0.5, 20000), ylim = c(0, 0.1), ylim.x = c(0, 1.2), f = 4)


par(mfrow = c(2, 4))

# pdf("All.pdf")

# pdf("AboveE.pdf")
aboveE(outer = FALSE, line = 1)

# pdf("Extra.pdf")
extra(outer = FALSE, line = 1)

# pdf("Above.pdf")
above(outer = FALSE, line = 1)

# pdf("Below.pdf")
below(outer = FALSE, line = 1)

# q()
