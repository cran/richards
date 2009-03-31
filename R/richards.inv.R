richards.inv <-
function (y, a = 0.1, d = 2.4, e = solveE(x50, b, ny), x50 = 100, 
    b = solveB(x = x50, a = a, d = d, ny = ny, x50 = x50, b4 = b4), 
    b4 = 1, ny = k - 1, k = 2) 
if (ny == 0) {
    f <- ifelse(min(a, d) < y & y < max(a, d), (y - d)/(a - d), 1)
    ifelse(min(a, d) < y & y < max(a, d), e * (-log(f))^(1/b), NA)
} else {
    ifelse(min(a, d) < y & y < max(a, d),
        e * (((((a - d)/(y - d))^(ny) - 1)/ny)^(1/b)), NA)
}

