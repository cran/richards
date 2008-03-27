`fpl` <-
function (x, a = 0.1, d = 2.4, e = 100, b = 1) 
d + (a - d)/(1 + (x/e)^b)
