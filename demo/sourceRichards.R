
sSrichards <- selfStart( ~ d + (a - d) / (1 + (2^ny-1)*(input/x50)^b)^(1/ny),
  function (mCall, data, LHS) initialRichards(mCall, data, LHS), c("a", "d", "b", "x50", "ny"))

sSrichardsE <- selfStart( ~ d + (a - d) / (1 + ny*(input/e)^b)^(1/ny),
  function (mCall, data, LHS) initialRichards(mCall, data, LHS, pnames = c("a", "d", "b", "e", "ny")), c("a", "d", "b", "e", "ny"))


#                           ~ d + (a - d) / (1 + (2^ny-1)*exp(-(xmid-input)/scal))^(1/ny),

s.richardsLog <- selfStart( ~ a + (d - a) / (1 + (2^ny-1)*exp((xmid-input)/scal))^(1/ny),
  function (mCall, data, LHS) initialRichardsLog(mCall, data, LHS), c("a", "d", "xmid", "scal", "ny"))

sSrichardsLog <- selfStart( ~ d + (a - d) / (1 + (2^ny-1)*exp((input-xmid)/scal))^(1/ny),
  function (mCall, data, LHS) initialRichardsLog(mCall, data, LHS), c("a", "d", "xmid", "scal", "ny"))


sSrichardsLogBM <- selfStart( ~ d + (a - d) / (1 + (2^(1/m)-1) * exp(b*(input-xmid)))^(m),
  function (mCall, data, LHS) initialRichardsLogB(mCall, data, LHS, 
            pnames = c("a", "d", "xmid", "b", "m")), c("a", "d", "xmid", "b", "m"))

sSnyFixedLogB <- selfStart( ~ d + (a - d) / (1 + (2^ny-1)*exp(b*(input-xmid)))^(1/ny),
  function (mCall, data, LHS) initialNyFixedLogB(mCall, data, LHS, ny = NULL), c("a", "d", "xmid", "b"))

# package.skeleton("SSrichards")
