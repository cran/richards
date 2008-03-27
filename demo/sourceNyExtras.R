
#          k        ny           factor         potens  
#  [1,]  -15       -16          -0.9999847     -1/16    
#  [2,]   -7        -8          -0.9960938     -1/8     
#  [3,]   -3        -4          -0.9375        -1/4     
#  [4,]   -1        -2          -3/4           -1/2     
#  [5,]    0        -1          -1/2           -1       
#  [6,]    1/5      -4/5        -0.4256508     -5/4     
#  [7,]    1/4      -3/4        -0.4053964     -4/3     
#  [8,]    1/3      -2/3        -0.3700395     -3/2     
#  [1,]    1/2      -1/2        -0.2928932     -2       
#  [2,]    2/3      -1/3        -0.2062995     -3       
#  [3,]    3/4      -1/4        -0.1591036     -4       
#  [4,]    4/5      -1/5        -0.1294494     -5       
#  [5,]    5/6      -1/6        -0.1091013     -6       
#  [6,]    1         0           0            Inf       
#  [7,]    7/6       1/6         0.122462       6       
#  [8,]    6/5       1/5         0.1486984      5       
#  [9,]    5/4       1/4         0.1892071      4       
# [10,]    4/3       1/3         0.2599211      3       
# [11,]    3/2       1/2         0.4142136      2       
#  [2,]    3         2           3              1/2     
#  [4,]    5         4          15              1/4     
#  [6,]    9         8         255              1/8     
#  [8,]   17        16       65535              1/16    
# [10,]   33        32  4294967295              1/32    
# [12,]   65        64           1.844674e+19   1/64    
# [14,]  129       128           3.402824e+38   1/128   

SS.d.nym1   <- selfStart( ~ d + (a - d) * (1 - 1/2 * (input/x50)^b),                        function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1  ), c("a", "d", "b", "x50"))
SS.d.nym1o3 <- selfStart( ~ d + (a - d) * (1 + (2^(-1/3)-1) * (input/x50)^b)^3,             function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1  ), c("a", "d", "b", "x50"))
SS.d.ny0    <- selfStart( ~ d + (a - d) * exp(- log(2) * (input / x50)^b),                  function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   0  ), c("a", "d", "b", "x50"))
SS.d.ny1    <- selfStart( ~ d + (a - d) / (1 + (input/x50)^b),                              function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1  ), c("a", "d", "b", "x50"))
SS.d.ny2    <- selfStart( ~ d + (a - d) / sqrt(1 + 3 * (input/x50)^b),                      function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   2  ), c("a", "d", "b", "x50"))
SS.d.ny3    <- selfStart( ~ d + (a - d) / (1 + 7 * (input/x50)^b)^(1/3),                    function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   3  ), c("a", "d", "b", "x50"))
SS.d.ny4    <- selfStart( ~ d + (a - d) / sqrt(sqrt((1 + 15 * (input/x50)^b))),             function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   4  ), c("a", "d", "b", "x50"))
SS.d.ny8    <- selfStart( ~ d + (a - d) / (1 + (255 - 1) * (input/x50)^b)^(1/8),            function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   8  ), c("a", "d", "b", "x50"))

sSny.m16   <- selfStart( ~ d + (a - d) / (1 + (2^(-16  ) - 1) * (input/x50)^b)^( -1/16  ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny = -16  ), c("a", "d", "b", "x50"))
sSny.m8    <- selfStart( ~ d + (a - d) / (1 + (2^( -8  ) - 1) * (input/x50)^b)^( -1/8   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -8  ), c("a", "d", "b", "x50"))
sSny.m4    <- selfStart( ~ d + (a - d) / (1 + (2^( -4  ) - 1) * (input/x50)^b)^( -1/4   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -4  ), c("a", "d", "b", "x50"))
sSny.m2    <- selfStart( ~ d + (a - d) / (1 +         - 3 / 4 * (input/x50)^b)^( -1/2   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -2  ), c("a", "d", "b", "x50"))
sSny.m1    <- selfStart( ~ d + (a - d) / (1 +         - 1 / 2 * (input/x50)^b)^( -1     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1  ), c("a", "d", "b", "x50"))
sSny.m4o5  <- selfStart( ~ d + (a - d) / (1 +  (2^(-4/5) - 1) * (input/x50)^b)^( -5/4   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -4/5), c("a", "d", "b", "x50"))
sSny.m3o4  <- selfStart( ~ d + (a - d) / (1 +  (2^(-3/4) - 1) * (input/x50)^b)^( -4/3   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -3/4), c("a", "d", "b", "x50"))
sSny.m2o3  <- selfStart( ~ d + (a - d) / (1 +  (2^(-2/3) - 1) * (input/x50)^b)^( -3/2   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -2/3), c("a", "d", "b", "x50"))
sSny.m1o2  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/2) - 1) * (input/x50)^b)^( -2     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/2), c("a", "d", "b", "x50"))
sSny.m1o3  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/3) - 1) * (input/x50)^b)^( -3     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/3), c("a", "d", "b", "x50"))
sSny.m1o4  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/4) - 1) * (input/x50)^b)^( -4     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/4), c("a", "d", "b", "x50"))
sSny.m1o5  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/5) - 1) * (input/x50)^b)^( -5     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/5), c("a", "d", "b", "x50"))
sSny.m1o6  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/6) - 1) * (input/x50)^b)^( -6     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/6), c("a", "d", "b", "x50"))
sSny.m1o7  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/7) - 1) * (input/x50)^b)^( -7     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/6), c("a", "d", "b", "x50"))
sSny.m1o8  <- selfStart( ~ d + (a - d) / (1 +  (2^(-1/8) - 1) * (input/x50)^b)^( -8     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/6), c("a", "d", "b", "x50"))
sSny.m1o16 <- selfStart( ~ d + (a - d) / (1 + (2^(-1/16) - 1) * (input/x50)^b)^( -16    ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/6), c("a", "d", "b", "x50"))
sSny.m1o32 <- selfStart( ~ d + (a - d) / (1 + (2^(-1/32) - 1) * (input/x50)^b)^( -32    ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  -1/6), c("a", "d", "b", "x50"))
# sSny.0   <- selfStart( ~ d + (a - d) / (1 +   (2^( 0 ) - 1) * (input/x50)^b)^( 1/0    ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   0  ), c("a", "d", "b", "x50"))
sSny.1o32  <- selfStart( ~ d + (a - d) / (1 + (2^( 1/32) - 1) * (input/x50)^b)^(  32     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/6), c("a", "d", "b", "x50"))
sSny.1o16  <- selfStart( ~ d + (a - d) / (1 + (2^( 1/16) - 1) * (input/x50)^b)^(  16    ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/6), c("a", "d", "b", "x50"))
sSny.1o8   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/8) - 1) * (input/x50)^b)^(  8     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/6), c("a", "d", "b", "x50"))
sSny.1o7   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/7) - 1) * (input/x50)^b)^(  7     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/6), c("a", "d", "b", "x50"))
sSny.1o6   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/6) - 1) * (input/x50)^b)^(  6     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/6), c("a", "d", "b", "x50"))
sSny.1o5   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/5) - 1) * (input/x50)^b)^(  5     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/5), c("a", "d", "b", "x50"))
sSny.1o4   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/4) - 1) * (input/x50)^b)^(  4     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/4), c("a", "d", "b", "x50"))
sSny.1o3   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/3) - 1) * (input/x50)^b)^(  3     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/3), c("a", "d", "b", "x50"))
sSny.1o2   <- selfStart( ~ d + (a - d) / (1 +  (2^( 1/2) - 1) * (input/x50)^b)^(  2     ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   1/2), c("a", "d", "b", "x50"))
sSny.1     <- selfStart( ~ d + (a - d) / (1 +               1 * (input/x50)^b)^(  1/1   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   2  ), c("a", "d", "b", "x50"))
sSny.2     <- selfStart( ~ d + (a - d) / (1 +               3 * (input/x50)^b)^(  1/2   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   2  ), c("a", "d", "b", "x50"))
sSny.3     <- selfStart( ~ d + (a - d) / (1 +               7 * (input/x50)^b)^(  1/3   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   2  ), c("a", "d", "b", "x50"))
sSny.4     <- selfStart( ~ d + (a - d) / (1 +              15 * (input/x50)^b)^(  1/4   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   4  ), c("a", "d", "b", "x50"))
sSny.8     <- selfStart( ~ d + (a - d) / (1 +             255 * (input/x50)^b)^(  1/8   ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =   8  ), c("a", "d", "b", "x50"))
sSny.16    <- selfStart( ~ d + (a - d) / (1 +           65535 * (input/x50)^b)^(  1/16  ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  16  ), c("a", "d", "b", "x50"))
sSny.32    <- selfStart( ~ d + (a - d) / (1 + (2^( 32  ) - 1) * (input/x50)^b)^(  1/32  ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  32  ), c("a", "d", "b", "x50"))
sSny.64    <- selfStart( ~ d + (a - d) / (1 + (2^( 64  ) - 1) * (input/x50)^b)^(  1/64  ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =  64  ), c("a", "d", "b", "x50"))
sSny.128   <- selfStart( ~ d + (a - d) / (1 + (2^(128  ) - 1) * (input/x50)^b)^(  1/128 ), function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny = 128  ), c("a", "d", "b", "x50"))


`SSny.m16` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^-16 - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -1/16
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
               initialNyFixed(mCall, data, LHS, ny = -16), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m8` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^-8 - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -1/8
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
               initialNyFixed(mCall, data, LHS, ny = -8), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m4` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^-4 - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -1/4
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
               initialNyFixed(mCall, data, LHS, ny = -4), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m2` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- -3/4
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -1/2
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
               initialNyFixed(mCall, data, LHS, ny = -2), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- -1 / 2
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -1
    .Z       <- 1 + .t2nym1 * .f502b
    .inv     <- .Z <= 0
    .Z2m     <- ifelse(.inv,     Inf, .Z^.m)
    .moZ     <- ifelse(.Z == 0,  NaN, .m / .Z)
    .rZ2m    <- ifelse(.inv,       0, 1 / .Z2m)
    .value   <- d + .aMd * .rZ2m
    .p5      <- ifelse(.inv,       0, .aMd * .t2nym1 * .rZ2m * .moZ * .f502b)
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "d", "b", "x50")))
    .grad[, "a"] <-      .rZ2m
    .grad[, "d"] <-  1 - .rZ2m
    .grad[, "b"] <-    - .p5 * log(.f50)
    .grad[, "x50"] <-    .p5 * b / x50
    attr(.value, "gradient") <- .grad
    .value
}, initial = function (mCall, data, LHS) 
               initialNyFixed(mCall, data, LHS, ny = -1), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m4o5` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-4/5) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -5/4
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
               initialNyFixed(mCall, data, LHS, ny = -4/5), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m3o4` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-3/4) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -4/3
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
               initialNyFixed(mCall, data, LHS, ny = -3/4), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m2o3` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-2/3) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -3/2
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
               initialNyFixed(mCall, data, LHS, ny = -2/3), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o2` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/2) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -2
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
               initialNyFixed(mCall, data, LHS, ny = -1/2), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o3` <-
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

`SSny.m1o4` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/4) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -4
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
               initialNyFixed(mCall, data, LHS, ny = -1/4), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o5` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/5) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -5
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
               initialNyFixed(mCall, data, LHS, ny = -1/5), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o6` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/6) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -6
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
               initialNyFixed(mCall, data, LHS, ny = -1/6), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o7` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/7) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -7
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
               initialNyFixed(mCall, data, LHS, ny = -1/7), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o8` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/8) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -8
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
               initialNyFixed(mCall, data, LHS, ny = -1/8), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o16` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/16) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -16
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
               initialNyFixed(mCall, data, LHS, ny = -1/16), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.m1o32` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(-1/32) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .m       <- -32
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
               initialNyFixed(mCall, data, LHS, ny = -1/32), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")


`SSny.1o32` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/32) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 32
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
initialNyFixed(mCall, data, LHS, ny = 1/32), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o16` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/16) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 16
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
initialNyFixed(mCall, data, LHS, ny = 1/16), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o8` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/8) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 8
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
initialNyFixed(mCall, data, LHS, ny = 1/8), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o7` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/7) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 7
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
initialNyFixed(mCall, data, LHS, ny = 1/7), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o6` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/6) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 6
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
initialNyFixed(mCall, data, LHS, ny = 1/6), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o5` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/5) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 5
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
initialNyFixed(mCall, data, LHS, ny = 1/5), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o4` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/4) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 4
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
initialNyFixed(mCall, data, LHS, ny = 1/4), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o3` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/3) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 3
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
initialNyFixed(mCall, data, LHS, ny = 1/3), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1o2` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(1/2) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 2
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
initialNyFixed(mCall, data, LHS, ny = 1/2), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.1` <-
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

`SSny.2` <-
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

`SSny.3` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 7 # 2^(3) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/3
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
initialNyFixed(mCall, data, LHS, ny = 3), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.4` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 15 # 2^(4) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/4
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
initialNyFixed(mCall, data, LHS, ny = 4), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.8` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(8) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/8
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
initialNyFixed(mCall, data, LHS, ny = 8), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.16` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(16) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/16
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
initialNyFixed(mCall, data, LHS, ny = 16), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.32` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(32) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/32
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
initialNyFixed(mCall, data, LHS, ny = 32), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.64` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(64) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/64
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
initialNyFixed(mCall, data, LHS, ny = 64), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")

`SSny.128` <-
structure(function (input, a, d, b, x50) 
{
    .aMd     <- a - d
    .t2nym1  <- 2^(128) - 1
    .f50     <- input / x50         # 'x50' (exp('xmid')) should not be 0!
    .f502b   <- .f50^b              # 'x' and 'x50' both positive (or negative)!
    .Z       <- 1 + .t2nym1 * .f502b
    .m       <- 1/128
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
initialNyFixed(mCall, data, LHS, ny = 128), 
 pnames = c("a", "d", "b", "x50"), class = "selfStart")
