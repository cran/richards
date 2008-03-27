`richards.deriv` <-
       function (x, a = 0.1, d = 2.4,
                 e = solveE(x50, b, ny), x50 = 100,
                 b = solveB(x = x50, a = a, d = d, ny = ny, 
                            x50 = x50, b4 = b4), b4 = 1,
                 ny = k - 1, k = 2)
if (ny == 0) {
    z <- x/e
    b * (d - a) * exp(-z^b) * (z^(b - 1))/e
} else {
    z <- x/e
    b * (d - a) * (z^(b - 1))/(e * ((1 + ny * z^b)^((1/ny) + 
        1)))
}
