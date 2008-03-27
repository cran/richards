`SSrichards` <-
structure(function (input, a, d, b, x50, ny) 
{
    .aMd     <- a - d
    .two2ny  <- 2^ny
    .t2nym1  <- .two2ny - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- 1 / ny              # 'ny' should not be 0!
    .Z       <- 1 + .t2nym1 * .f502b
    .inv     <- .Z <= 0
    .log.Z   <- log(ifelse(.inv, NaN, .Z))
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p1      <- .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .p1 * .t2nym1 * .moZ * .f502b)
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a",  "d", "b", "x50", "ny")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "b"]     <-   - .p5 * log(.f50)
    # -((a-d)*(2^ny-1)*log(x/x50)*(x/x50)^b*Z^(-1/ny-1))/ny
    # = - .aMd * .t2nym1 * log(.f50) * .f502b * .rZ2m * .moZ
    .grad[, "x50"]   <-     .p5 * b / x50
    # (b*(a-d)*(2^ny-1)*(x/x50)^b*Z^(-1/ny-1))/(ny*x50)
    # = b * .aMd * .t2nym1 * .f502b * .rZ2m * .moZ / x50
    .grad[, "ny"]    <- ifelse(.inv,     0,
            .p1 * (.log.Z * .m^2 - log(2) * .two2ny * .f502b * .moZ))
    # ((a-d)*(log(Z)/ny^2-(log(2)*2^ny*(x/x50)^b)/(ny*Z)))/Z^(1/ny)
    # = .aMd * (.log.Z * .m^2 - log(2) * .two2ny * .f502b * .moZ)) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichards(mCall, data, LHS),
 pnames = c("a", "d", "b", "x50", "ny"), class = "selfStart")
