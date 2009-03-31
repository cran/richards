`SSrichardsLogBG` <-
structure(function (input, a, d, xmid, b, g) 
{
    .aMd     <- a - d
    .two2ny  <- 2^(1/g)             # 'g' should not be 0!
    .t2nym1  <- .two2ny - 1
    .dmid    <- input - xmid
    .edmids  <- exp(b * .dmid)
    .Z       <- 1 + .t2nym1 * .edmids
    .inv     <- .Z <= 0
    .log.Z   <- log(ifelse(.inv, NaN, .Z))
    .Z2m     <- ifelse(.inv,     Inf, .Z^g)
    .rZ      <- ifelse(.Z == 0,  NaN, 1 / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p1      <- .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .p1 * .t2nym1 * g * .rZ * .edmids)
    .value   <- d + .aMd/.Z2m
    .grad    <- array(0, c(length(.value), 5L), 
                      list(NULL, c("a", "d", "xmid", "b", "g")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "xmid"]  <-     .p5 * b
    # b*(a-d)*g*(2^(1/g)-1)*%e^(b*(input-xmid))*Z^(-g-1)
    # b * .aMd * g * .t2nym1 * .edmids * .rZ2m * .rZ
    .grad[, "b"]     <-   - .p5 * .dmid
    # -(a-d)*g*(2^(1/g)-1)*(input-xmid)*%e^(b*(input-xmid))*Z^(-g-1)
    # - .aMd * g * .t2nym1 * .dmid * .edmids * .rZ2m * .rZ
    .grad[, "g"]     <- ifelse(.inv,     0,
            .p1 * ( .two2ny * log(2) * .edmids * .rZ / g - .log.Z))
    # ((a-d)*((log(2)*2^(1/g)*%e^(b*(input-xmid)))/(g*Z)-log(Z)))/Z^g
    # .aMd * ( log(2) * .two2ny * .edmids * (1/g) * .rZ - .log.Z) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichardsLogB(mCall, data, LHS, 
                                   pnames = c("a", "d", "xmid", "b", "g")), 
 pnames = c("a", "d", "xmid", "b", "g"), class = "selfStart")
