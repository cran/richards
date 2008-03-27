`SSnym1o3` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/3) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -3
    .Z       <- 1 + .t2nym1 * .f502b
    .inv     <- .Z <= 0
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .aMd * .t2nym1 * .rZ2m * .moZ * .f502b)
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "b"]     <-   - .p5 * log(.f50)
    .grad[, "x50"]   <-     .p5 * b / x50
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialNyFixed(mCall, data, LHS, ny = -1/3), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")
