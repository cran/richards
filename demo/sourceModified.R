`sSnyFixedLog.Cond` <-
structure(function (ny, input, a, d, xmid, scal) 
{
    .expr1  <- a - d
    .expr3  <- 2^ny - 1
    .expr4  <- input - xmid
    .expr6  <- exp(.expr4/scal)
    .expr8  <- 1 + .expr3 * .expr6
    .expr9  <- 1/ny
    .inv    <- .expr8 <= 0
    .expr10 <- ifelse(.inv,   Inf, .expr8^.expr9)
    .expr13 <- ifelse(.inv,     0, 1 / .expr10)
    .expr16 <- ifelse(.inv, - Inf, .expr8^(.expr9 - 1))
    .expr23 <- ifelse(.inv,   Inf, .expr10^2)
    .value  <- ifelse(.inv, d + 0, d + .expr1 / .expr10)
    .expr31 <- ifelse(.inv,     0,
                      .expr1 * .expr16 * .expr9 * .expr3 * .expr6 / .expr23)
    .grad <- array(0, c(length(.value), 4L),
                   list(NULL, c("a", "d", "xmid", "scal")))
    .grad[, "a"]    <- .expr13
    .grad[, "d"]    <- 1 - .expr13
    .grad[, "xmid"] <- .expr31 / scal
    .grad[, "scal"] <- .expr31 * .expr4 / scal^2
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialNyFixedLog(mCall, data, LHS, ny = NULL),
 pnames = c("a", "d", "xmid", "scal"), class = "selfStart")

`sSnyFixed.Cond` <-
structure(function (ny, input, a, d, b, x50) 
{
    .expr1  <- a - d
    .expr3  <- 2^ny - 1
    .expr4  <- input/x50
    .expr5  <- .expr4^b
    .expr7  <- 1 + .expr3 * .expr5
    .expr8  <- 1/ny
    .inv    <- .expr7 <= 0
    .expr9  <- ifelse(.inv,   Inf, .expr7^.expr8)
    .expr12 <- ifelse(.inv,     0, 1 / .expr9)
    .expr15 <- ifelse(.inv, - Inf, .expr7^(.expr8 - 1))
    .expr22 <- ifelse(.inv,   Inf, .expr9^2)
    .value  <- ifelse(.inv, d + 0, d + .expr1 / .expr9)
    .expr31 <- ifelse(.inv,     0,
                      .expr1 * .expr15 * .expr8 * .expr3 / .expr22)
    .grad   <- array(0, c(length(.value), 4L), 
                     list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"]   <- .expr12
    .grad[, "d"]   <- 1 - .expr12
    .grad[, "b"]   <-   - .expr31 * .expr5 * log(.expr4)
    .grad[, "x50"] <-     .expr31 * .expr4^(b - 1) * b * (input/x50^2)
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialNyFixed(mCall, data, LHS, ny = NULL),
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`sSrichardsLog.Cond` <-
structure(function (input, a, d, xmid, scal, ny) 
{
    .expr1   <- a - d
    .expr2   <- 2^ny
    .expr3   <- .expr2 - 1
    .expr4   <- input - xmid
    .expr6   <- exp(.expr4/scal)
    .expr8   <- 1 + .expr3 * .expr6
    .expr9   <- 1/ny
    .inv     <- .expr8 <= 0
    .expr10  <- ifelse(.inv,   Inf, .expr8^.expr9)
    .expr13  <- ifelse(.inv,     0, 1 / .expr10)
    .expr16  <- ifelse(.inv, - Inf, .expr8^(.expr9 - 1))
    .expr23  <- ifelse(.inv,   Inf, .expr10^2)
    .value   <- ifelse(.inv, d + 0, d + .expr1  / .expr10)
    .expr123 <- ifelse(.inv,     0, .expr1  / .expr23)
    .expr169 <- ifelse(.inv,     0, .expr16 * .expr9)
    .expr36  <- .expr123 * .expr169 * .expr3 * .expr6
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a", "d", "xmid", "scal", "ny")))
    .grad[, "a"]    <- .expr13
    .grad[, "d"]    <- 1 - .expr13
    .grad[, "xmid"] <- .expr36 / scal
    .grad[, "scal"] <- .expr36 * .expr4/scal^2
    .log.expr8      <- log(ifelse(.inv, 1, .expr8))
    .grad[, "ny"]   <- ifelse(.inv,     0,
                       - .expr123 * (.expr169 * .expr2 * log(2) * .expr6
                                    - .expr10 * .log.expr8 / ny^2))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichardsLog(mCall, data, LHS),
 pnames = c("a", "d", "xmid", "scal", "ny"), class = "selfStart")

`sSrichards.Cond` <-
structure(function (input, a, d, b, x50, ny) 
{
    .expr1   <- a - d
    .expr2   <- 2^ny
    .expr3   <- .expr2 - 1
    .expr4   <- input/x50
    .expr5   <- .expr4^b
    .expr7   <- 1 + .expr3 * .expr5
    .expr8   <- 1/ny
    .inv     <- .expr7 <= 0
    .expr9   <- ifelse(.inv,   Inf, .expr7^.expr8)
    .expr12  <- ifelse(.inv,     0, 1 / .expr9)
    .expr15  <- ifelse(.inv, - Inf, .expr7^(.expr8 - 1))
    .expr22  <- ifelse(.inv,   Inf, .expr9^2)
    .value   <- ifelse(.inv, d + 0, d + .expr1 / .expr9)
    .expr122 <- ifelse(.inv,     0, .expr1 / .expr22)
    .expr158 <- ifelse(.inv,     0, .expr15 * .expr8)
    .expr35  <- .expr122 * .expr158 * .expr3
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a",  "d", "b", "x50", "ny")))
    .grad[, "a"]   <- .expr12
    .grad[, "d"]   <- 1 - .expr12
    .grad[, "b"]   <- - .expr35 * .expr5 * log(.expr4)
    .grad[, "x50"] <-   .expr35 * .expr4^(b - 1) * b * (input/x50^2)
    .log.expr7     <- log(ifelse(.inv, 1, .expr7))
    .grad[, "ny"]  <- ifelse(.inv,     0,
                      - .expr122 * (.expr158 * .expr2 * log(2) * .expr5
                                   - .expr9 * .log.expr7 / ny^2))
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichards(mCall, data, LHS),
 pnames = c("a", "d", "b", "x50", "ny"), class = "selfStart")

`sSrichardsE.Cond` <-
structure(function (input, a, d, b, e, ny) 
{
    .aMd     <- a - d
    .expr2   <- input/e
    .two2ny   <- .expr2^b
    .Z       <- 1 + ny * .two2ny
    .m       <- 1/ny
    .inv     <- .Z <= 0
    .Z2m     <- ifelse(.inv,   Inf, .Z^.m)
    .rZ2m    <- ifelse(.inv,     0, 1/.Z2m)
    .Z2mm1   <- ifelse(.inv, - Inf, .Z^(.m - 1))
    .rZ2ms   <- ifelse(.inv,     0, .rZ2m^2)
    .value   <- ifelse(.inv,     d, d + .aMd/.Z2m)
    .p1      <- ifelse(.inv,     0, .aMd * .rZ2ms)
    .p2m     <- ifelse(.inv,     0, .Z2mm1 * .m)
    .p5      <- .p1 * .p2m * ny
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a", "d", "b", "e", "ny")))
    .grad[, "a"]  <- .rZ2m
    .grad[, "d"]  <- 1 - .rZ2m
    .grad[, "b"]  <- - .p5 * .two2ny * log(.expr2)
    .grad[, "e"]  <-   .p5 * .expr2^(b - 1) * b * (input/e^2)
    .log.Z    <- log(ifelse(.inv, 1, .Z))
    .grad[, "ny"] <-   .p1 * (.Z2m * .log.Z * (1/ny^2) - .p2m * .two2ny)
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
initialRichards(mCall, data, LHS,
 pnames = c("a", "d", "b", "e", "ny")),
 pnames = c("a", "d", "b", "e", "ny"), class = "selfStart")

`sS.nyFixedLog` <-
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

`sS.nyFixed` <-
structure(function (ny, input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^ny - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- 1 / ny              # 'ny' should not be 0!
    .Z       <- 1 + .t2nym1 * .f502b
    .inv     <- .Z <= 0
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .aMd * .t2nym1 * .rZ2m * .moZ * .f502b)
    .grad    <- array(0, c(length(.value), 4L), 
                      list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "b"]     <-   - .p5 * log(.f50)
    .grad[, "x50"]   <-     .p5 * b / x50
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialNyFixed(mCall, data, LHS, ny = NULL),
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`sS.richardsLog` <-
structure(function (input, a, d, xmid, scal, ny) 
{
    .aMd     <- a - d
    .two2ny  <- 2^ny
    .t2nym1  <- .two2ny - 1
    .dmid    <- input - xmid
    .rscal   <- 1 / scal
    .edmids  <- exp(.dmid * .rscal) # 'scal' should not be 0!
    .m       <- 1 / ny              # 'ny' should not be 0!
    .Z       <- 1 + .t2nym1 * .edmids
    .inv     <- .Z <= 0
    .log.Z   <- log(ifelse(.inv, NaN, .Z))
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p1      <- .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .p1 * .t2nym1 * .moZ * .edmids)
    .grad    <- array(0, c(length(.value), 5L),
                      list(NULL, c("a", "d", "xmid", "scal", "ny")))
    .grad[, "a"]     <-     .rZ2m
    .grad[, "d"]     <- 1 - .rZ2m
    .grad[, "xmid"]  <-     .p5 * .rscal
    # ((a-d)*(2^ny-1)*%e^((input-xmid)/scal)*Z^(-1/ny-1))/(ny*scal)
    # = .aMd * .t2nym1 * .edmids * .rZ2m * .moZ * .rscal
    .grad[, "scal"]  <-     .p5 * .dmid * .rscal^2
    # ((a-d)*(2^ny-1)*(input-xmid)*%e^((input-xmid)/scal)*Z^(-1/ny-1))/(ny*scal^2)
    # = .aMd * .t2nym1 * .dmid * .edmids * .rZ2m * .moZ * .rscal^2
    .grad[, "ny"]    <- ifelse(.inv,     0,
            .p1 * (.log.Z * .m^2 - log(2) * .two2ny * .edmids * .moZ))
    # ((a-d)*(log(Z)/ny^2-(log(2)*2^ny*%e^((input-xmid)/scal))/(ny*Z)))/Z^(1/ny)
    # = .aMd * (.log.Z * .m^2 - log(2) * .two2ny * .edmids * .moZ)) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichardsLog(mCall, data, LHS),
 pnames = c("a", "d", "xmid", "scal", "ny"), class = "selfStart")

`sS.richardsLogBG` <-
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
            .p1 * ( .two2ny * log(2) * .edmids * .rZ / g - log(.Z)))
    # ((a-d)*((log(2)*2^(1/g)*%e^(b*(input-xmid)))/(g*Z)-log(Z)))/Z^g
    # .aMd * ( log(2) * .two2ny * .edmids * (1/g) * .rZ - log(.Z)) * .rZ2m
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialRichardsLogB(mCall, data, LHS, 
                                   pnames = c("a", "d", "xmid", "b", "g")), 
 pnames = c("a", "d", "xmid", "b", "g"), class = "selfStart")

`sS.richards` <-
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

`sS.richardsE` <-
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

