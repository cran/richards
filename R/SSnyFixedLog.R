`SSnyFixedLog` <-
structure(function (ny, input, a, d, xmid, scal) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^ny - 1
    .dmid    <- input - xmid
    .rscal   <- 1 / scal
    .edmids  <- exp(.dmid * .rscal) # 'scal' should not be 0!
    .m       <- 1 / ny              # 'ny' should not be 0!
    .Z       <- 1 + .t2nym1 * .edmids
    .inv     <- .Z <= 0
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .aMd * .t2nym1 * .rZ2m * .moZ * .edmids)
    .grad    <- array(0, c(length(.value), 4L),
                      list(NULL, c("a", "d", "xmid", "scal")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "xmid"]  <-     .p5 * .rscal
    .grad[, "scal"]  <-     .p5 * .dmid * .rscal^2
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialNyFixedLog(mCall, data, LHS, ny = NULL),
 pnames = c("a", "d", "xmid", "scal"), class = "selfStart")
