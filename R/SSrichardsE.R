`SSrichardsE` <-
structure(function (input, a, d, b, e, ny) 
{
    .aMd     <- a - d
    .fe      <- input / e           # 'e' (~ exp('xmid')) should not be 0!
    .fe2b    <- .fe^b               # 'x' and 'e' both positive (or negative)!
    .m       <- 1 / ny              # 'ny' should not be 0!
    .Z       <- 1 + ny * .fe2b
    .inv     <- .Z <= 0
    .log.Z   <- log(ifelse(.inv, NaN, .Z))
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .rZ      <- ifelse(.Z == 0,  NaN, 1 / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p1      <- .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .p1 * .rZ * .fe2b)
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a", "d", "b", "e", "ny")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "b"]     <-   - .p5 * log(.fe)
    # -(a-d)*(x/e)^b*log(x/e)*Z^(-1/ny-1)
    # = - .aMd * .fe2b * log(.fe) * .rZ2m * .rZ 
    .grad[, "e"]     <-     .p5 * b / e
    # (b*(a-d)*(x/e)^b*Z^(-1/ny-1))/e
    # = b * .aMd * .fe2b * .rZ2m * .rZ / e
    .grad[, "ny"]    <- ifelse(.inv,     0, 
                              .p1 * (.log.Z * .m^2 - .rZ * .m * .fe2b))
    # ((a-d)*(log(Z)/ny^2-(x/e)^b/(ny*Z)))/Z^(1/ny)
    # = .aMd * (.log.Z * .m^2 - .fe2b * .rZ * .m) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialRichards(mCall, data, LHS,
 pnames = c("a", "d", "b", "e", "ny")),
 pnames = c("a", "d", "b", "e", "ny"), class = "selfStart")
