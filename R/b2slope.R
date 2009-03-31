`b2slope` <-
function (ny = X[, "ny"], b = X[, "b"], a = X[, "a"], d = X[, "d"],
          X = NULL, fold = 2, delta = 1e-09) 
    b * (d - a) * ny2f(ny, delta = delta) * log(fold)
