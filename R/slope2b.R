`slope2b` <-
function (ny = X[, "ny"], slope = X[, "slope"], a = X[, "a"], d = X[, "d"], 
          X = NULL, fold = 2, delta = 1e-09) 
    slope/((d - a) * ny2f(ny, delta = delta) * log(fold))
