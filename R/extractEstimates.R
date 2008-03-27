`extractEstimates` <-
function (fits, n = 5) 
{
    estimates <- lapply(fits, function(i) if (is.null(i) | class(i) == 
        "try-error") 
        rep(NA, n)
    else summary(i)$parameters[, "Estimate"])
    result <- matrix(unlist(estimates), ncol = n, byrow = TRUE)
    dimnames(result) <- list(names(fits), unlist(lapply(estimates, 
        names))[1:n])
    return(result)
}
