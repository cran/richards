`solveE` <-
function (x50, b, ny = k - 1, k = 2) 
x50/ifelse(ny == 0, log(2), abs((2^(ny) - 1)/ny))^(1/b)
