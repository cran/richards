`SSrichardsLog` <-
structure(function (input, a, d, xmid, scal, ny) 
{
    .aMd     <- a - d
    .two2ny  <- 2^ny
    .t2nym1  <- .two2ny - 1
    .dmid    <- input - xmid
    .rscal   <- 1 / scal
    .edmids  <- exp(.dmid * .rscal) # 'scal' should not be 0!
    .m       <- 1 / ny              # 'ny' should not be 0!
    .Z       <- 1 + .t2nym1 * .edmids
    .inv     <- .Z <= 0
    .log.Z   <- log(ifelse(.inv, NaN, .Z))
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p1      <- .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .p1 * .t2nym1 * .moZ * .edmids)
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a", "d", "xmid", "scal", "ny")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "xmid"]  <-     .p5 * .rscal
    # ((a-d)*(2^ny-1)*%e^((input-xmid)/scal)*Z^(-1/ny-1))/(ny*scal)
    # = .aMd * .t2nym1 * .edmids * .rZ2m * .moZ * .rscal
    .grad[, "scal"]  <-     .p5 * .dmid * .rscal^2
    # ((a-d)*(2^ny-1)*(input-xmid)*%e^((input-xmid)/scal)*Z^(-1/ny-1))/(ny*scal^2)
    # = .aMd * .t2nym1 * .dmid * .edmids * .rZ2m * .moZ * .rscal^2
    .grad[, "ny"]    <- ifelse(.inv,     0,
            .p1 * (.log.Z * .m^2 - log(2) * .two2ny * .edmids * .moZ))
    # ((a-d)*(log(Z)/ny^2-(log(2)*2^ny*%e^((input-xmid)/scal))/(ny*Z)))/Z^(1/ny)
    # = .aMd * (.log.Z * .m^2 - log(2) * .two2ny * .edmids * .moZ)) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichardsLog(mCall, data, LHS),
 pnames = c("a", "d", "xmid", "scal", "ny"), class = "selfStart")
