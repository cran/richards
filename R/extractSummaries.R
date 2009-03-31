`extractSummaries` <-
function (fits) 
{
    invalidFit <- function(fit)
        is.null(fit) | any(class(fit) == c("try-error", "data-error"))
    result <- lapply(fits, function(i) if (invalidFit(i)) 
        NULL
    else summary(i))
    return(result)
}
