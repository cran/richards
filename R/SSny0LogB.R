`SSny0LogB` <-
structure(function (input, a, d, xmid, b) 
{
    .aMd <- a - d
    .log2 <- log(2)
    .dmid <- input - xmid
    .edmids <- exp(b * .dmid)
    .ee <- exp(- .log2 * .edmids)
    .value <- d + .aMd * .ee
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "xmid", "b")))
    .grad[, "a"] <- .ee
    .grad[, "d"] <- 1 - .ee
    .grad[, "xmid"] <- .aMd * (.ee * (.log2* (.edmids * b)))
    .grad[, "b"] <- -(.aMd * (.ee * (.log2* (.edmids * .dmid))))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialNyFixedLogB(mCall, data, LHS, ny = 0), 
 pnames = c("a", "d", "xmid", "b"), class = "selfStart")
