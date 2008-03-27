
       sSny1 <- selfStart( ~ d + (a - d) / (1 + (input/x50)^b),                            
  function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =     1), c("a", "d", "b", "x50"))
       sSny2 <- selfStart( ~ d + (a - d) / sqrt(1 + 3 * (input/x50)^b),                    
  function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =     2), c("a", "d", "b", "x50"))
       sSny3 <- selfStart( ~ d + (a - d) / (1 + 7 * (input/x50)^b)^(1/3),                  
  function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =     3), c("a", "d", "b", "x50"))
       sSny4 <- selfStart( ~ d + (a - d) / sqrt(sqrt((1 + 15 * (input/x50)^b))),           
  function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =     4), c("a", "d", "b", "x50"))

   sSnyFixed <- selfStart( ~ d + (a - d) / (1 +        (2^ny-1) * (input/x50)^b)^(  1/ny  ),
  function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny = NULL ), c("a", "d", "b", "x50"))

s.nyFixedLog <- selfStart( ~ a + (d - a) / (1 + (2^ny-1)*exp((xmid-input)/scal))^(1/ny),
  function (mCall, data, LHS) initialNyFixedLog(mCall, data, LHS, ny = NULL), c("a", "d", "xmid", "scal"))

sSnyFixedLog <- selfStart( ~ d + (a - d) / (1 + (2^ny-1)*exp((input-xmid)/scal))^(1/ny),
  function (mCall, data, LHS) initialNyFixedLog(mCall, data, LHS, ny = NULL), c("a", "d", "xmid", "scal"))

sSnyFixedLogB <- selfStart( ~ d + (a - d) / (1 + (2^ny-1)*exp(b * (input-xmid)))^(1/ny),
  function (mCall, data, LHS) initialNyFixedLogB(mCall, data, LHS, ny = NULL), c("a", "d", "xmid", "b"))


 # (input / x50)^b ->

 # (exp(input) / exp(x50))^b =
 # exp(log((exp(input) / exp(x50))^b) =
 # exp(b * (input - xmid)) =
 # exp(b * input) * exp( - b * xmid) =
 # exp(log(b3) * input) * exp( - b * xmid) =
 # b3^input * exp(xmid) ^ (-b) =
 # b3^input * (log(2) * e ^ (-b)

 # b3 = exp(b)
 # x50 = exp(xmid) = e * log(2) ^ (1/b)

     sSny0 <- selfStart( ~ d + (a - d) * exp(- log(2) * (input / x50)^b),
   function (mCall, data, LHS) initialNyFixed(mCall, data, LHS, ny =     0), c("a", "d", "b", "x50"))

  sSny0Log <- selfStart( ~ d + (a - d) * exp(- log(2) *  exp((input - xmid) / scal)),
   function (mCall, data, LHS) initialNyFixedLog(mCall, data, LHS, ny = 0), c("a", "d", "xmid", "scal"))

 sSny0LogB <- selfStart( ~ d + (a - d) * exp(- log(2) *  exp(b * (input - xmid))),
   function (mCall, data, LHS) initialNyFixedLogB(mCall, data, LHS, ny = 0), c("a", "d", "xmid", "b"))

