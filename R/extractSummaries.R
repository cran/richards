`extractSummaries` <-
function (fits) 
{
    result <- lapply(fits, function(i) if (is.null(i) | class(i) == 
        "try-error") 
        NULL
    else summary(i))
    return(result)
}
