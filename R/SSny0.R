`SSny0` <-
structure(function (input, a, d, b, x50) 
{
    .expr1 <- a - d
    .expr2 <- log(2)
    .expr4 <- input/x50
    .expr5 <- .expr4^b
    .expr7 <- exp(-.expr2 * .expr5)
    .value <- d + .expr1 * .expr7
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", 
        "d", "b", "x50")))
    .grad[, "a"] <- .expr7
    .grad[, "d"] <- 1 - .expr7
    .grad[, "b"] <- -(.expr1 * (.expr7 * (.expr2 * (.expr5 * 
        log(.expr4)))))
    .grad[, "x50"] <- .expr1 * (.expr7 * (.expr2 * (.expr4^(b - 
        1) * (b * (input/x50^2)))))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialNyFixed(mCall, data, LHS, ny = 0), pnames = c("a", "d", 
"b", "x50"), class = "selfStart")


`SSny0` <-
structure(function (input, a, d, b, x50) 
{
    .aMd <- a - d
    .log2 <- log(2)
    .f50 <- input/x50
    .f502b <- .f50^b
    .ee <- exp(- .log2 * .f502b)
    .value <- d + .aMd * .ee
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"] <- .ee
    .grad[, "d"] <- 1 - .ee
    .grad[, "b"] <- -(.aMd * (.ee * (.log2 * (.f502b * log(.f50)))))
    .grad[, "x50"] <- .aMd * (.ee * (.log2 * (.f50^(b - 1) * (b * (input/x50^2)))))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialNyFixed(mCall, data, LHS, ny = 0), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")
