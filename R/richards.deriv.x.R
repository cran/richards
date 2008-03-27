`richards.deriv.x` <-
       function (x, a = 0.1, d = 2.4,
                 e = solveE(x50, b, ny), x50 = 100,
                 b = solveB(x = x50, a = a, d = d, ny = ny, 
                            x50 = x50, b4 = b4), b4 = 1,
                 ny = k - 1, k = 2)
richards.deriv(x, a = a, d = d, e = e, b = b, ny = ny) * x
