`SSny1` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .f502b
    .rZ2m    <- 1 / .Z # 1
    .Zs      <- .Z^2
    .value   <- d + .aMd / .Z # 1
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"] <-     .rZ2m
    .grad[, "d"] <- 1 - .rZ2m
    .grad[, "b"] <-   - .aMd / .Zs * .f502b * log(.f50) # 1
    .grad[, "x50"] <-   .aMd / .Zs * .f502b * b /x50    # 1
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialNyFixed(mCall, data, LHS, ny = 1), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")
