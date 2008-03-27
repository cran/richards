`fitList` <-
function (list, FUN,
          applyFUN2X = function(X, FUN, ...) FUN(X, ...),
          n = 4, ...) 
{
    fits <- lapply(list, function(X) applyFUN2X(X, FUN, n = n, 
        ...))
    estimates <- extractEstimates(fits, n = n)
    print(estimates)
    return(list(FUN = FUN, fits = fits, estimates = estimates))
}
