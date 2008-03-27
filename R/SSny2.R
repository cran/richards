`SSny2` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 3 # 2^(2) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/2
    .Z2m     <- .Z^.m
    .rZ2m    <- 1 / .Z2m
    .value   <- d + .aMd * .rZ2m
    .p5      <- .aMd * .t2nym1 * .rZ2m * .m / .Z * .f502b
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"] <-     .rZ2m
    .grad[, "d"] <- 1 - .rZ2m
    .grad[, "b"] <-   - .p5 * log(.f50)
    .grad[, "x50"] <-   .p5 * b / x50
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialNyFixed(mCall, data, LHS, ny = 2), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")
