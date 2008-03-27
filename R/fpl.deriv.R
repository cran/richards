`fpl.deriv` <-
function (x, a = 0.1, d = 2.4, e = 100, b = 1) 
if (b == 1) -(a - d) * (1/e)/(1 + (x/e)^b)^2 else -(a - d) * 
    b * (x/e)^(b - 1)/(1 + (x/e)^b)^2/e
