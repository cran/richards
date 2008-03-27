`SSrichardsLogBM` <-
structure(function (input, a, d, xmid, b, m) 
{
    .aMd     <- a - d
    .two2ny  <- 2^(1/m)             # 'm' should not be 0!
    .t2nym1  <- .two2ny - 1
    .dmid    <- input - xmid
    .edmids  <- exp(b * .dmid)
    .Z       <- 1 + .t2nym1 * .edmids
    .inv     <- .Z <= 0
    .log.Z   <- log(ifelse(.inv, NaN, .Z))
    .Z2m     <- ifelse(.inv,     Inf, .Z^m)
    .rZ      <- ifelse(.Z == 0,  NaN, 1 / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p1      <- .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .p1 * .t2nym1 * m * .rZ * .edmids)
    .value   <- d + .aMd/.Z2m
    .grad    <- array(0, c(length(.value), 5L), 
                      list(NULL, c("a", "d", "xmid", "b", "m")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "xmid"]  <-     .p5 * b
    # b*(a-d)*m*(2^(1/m)-1)*%e^(b*(input-xmid))*Z^(-m-1)
    # b * .aMd * m * .t2nym1 * .edmids * .rZ2m * .rZ
    .grad[, "b"]     <-   - .p5 * .dmid
    # -(a-d)*m*(2^(1/m)-1)*(input-xmid)*%e^(b*(input-xmid))*Z^(-m-1)
    # - .aMd * m * .t2nym1 * .dmid * .edmids * .rZ2m * .rZ
    .grad[, "m"]     <- ifelse(.inv,     0,
            .p1 * ( .two2ny * log(2) * .edmids * .rZ / m - .log.Z))
    # ((a-d)*((log(2)*2^(1/m)*%e^(b*(input-xmid)))/(m*Z)-log(Z)))/Z^m
    # .aMd * ( log(2) * .two2ny * .edmids * (1/m) * .rZ - .log.Z) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichardsLogB(mCall, data, LHS, 
                                   pnames = c("a", "d", "xmid", "b", "m")), 
 pnames = c("a", "d", "xmid", "b", "m"), class = "selfStart")
