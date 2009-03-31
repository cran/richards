`extractEstimates` <-
function (fits, n = 5) 
{
    invalidFit <- function(fit)
        is.null(fit) | any(class(fit) == c("try-error", "data-error"))
    estimates <- lapply(fits, function(i)
        if (invalidFit(i))
        rep(NA, n)
    else summary(i)$parameters[, "Estimate"])
    result <- matrix(unlist(estimates), ncol = n, byrow = TRUE)
    dimnames(result) <- list(names(fits), unlist(lapply(estimates, 
        names))[1:n])
    return(result)
}
