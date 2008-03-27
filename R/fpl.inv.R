`fpl.inv` <-
function (y, a = 0.1, d = 2.4, e = 100, b = 1) 
ifelse(min(a, d) < y & y < max(a, d), e * (((a - y)/(y - d))^(1/b)), NA)
