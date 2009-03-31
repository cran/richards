`x2c` <-
function (x50 = X[, "x50"], ny = X[, "ny"], b = X[, "b"], X = NULL, 
          fold = 2, x50f = 1, delta = 1e-09) 
{
    ny2logz50b <- function(ny, delta = 1e-09) ifelse(abs(ny) < 
        delta, log(log(2)), log((exp(log(2) * ny) - 1)/ny))
    logz50b <- ny2logz50b(ny, delta)
    E <- x50f * x50/exp(logz50b/b)
    h <- ifelse(abs(ny) < delta, 1, exp(log(abs(ny))/b))
    C <- x50f * x50/exp(logz50b/b)/h
    logfC <- log(x50f * x50/exp(logz50b/b)/h)/log(fold)
    return(cbind(logz50b, E, h, C, logfC))
}
