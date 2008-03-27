`SSny0Log` <-
structure(function (input, a, d, xmid, scal) 
{
    .aMd <- a - d
    .log2 <- log(2)
    .dmid <- input - xmid
    .edmids <- exp(.dmid/scal)
    .ee <- exp(- .log2* .edmids)
    .value <- d + .aMd * .ee
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "xmid", "scal")))
    .grad[, "a"] <- .ee
    .grad[, "d"] <- 1 - .ee
    .grad[, "xmid"] <- .aMd * (.ee * (.log2* (.edmids * (1/scal))))
    .grad[, "scal"] <- .aMd * (.ee * (.log2* (.edmids * (.dmid/scal^2))))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialNyFixedLog(mCall, data, LHS, ny = 0), 
 pnames = c("a", "d", "xmid", "scal"), class = "selfStart")
