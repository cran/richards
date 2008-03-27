`solveB` <-
function (x = x50, a = 0.1, d = 2.4, x50 = 100, b4 = 1, ny = k - 1, k = 2) 
{
    iterateB <- function(x, a = 0.1, d = 2.4, 
                         e = solveE(x50, b, ny), x50 = 100, b = 1, 
                         ny = 1 - k, k = 2, fd = 1, z = x/e) {
        if (ny == 0) {
            r <- fd * x/((d - a) * z^b * exp(-z^b))
        }
        else {
            r <- fd/((d - a) * (z^(b - 1))/(e * ((1 + ny * z^b)^((1/ny) + 1))))
        }
        return(r)
    }
    fd <- fpl.deriv(x, a = a, d = d, e = x50, b = b4)
    beta <- b4
    for (i in 1:5) beta <- iterateB(x, a = a, d = d, x50 = x50, b = beta, 
                                    ny = ny, fd = fd)
    return(beta)
}
