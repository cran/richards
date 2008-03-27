`fitDoubleList` <-
function (list, FUN, 
          applyFUN2X = function(X, FUN, ...) FUN(X, ...),
          n = 4, ...) 
{
    fits <- lapply(list, function(i) lapply(i, function(X) applyFUN2X(X, 
        FUN, n = n, ...)))
    estimates <- lapply(fits, function(i) extractEstimates(i, 
        n = n))
    return(list(FUN = FUN, fits = fits, estimates = estimates))
}
