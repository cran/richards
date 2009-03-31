`ny2f` <-
function (ny, delta = 1e-09) 
{
    ny2logz50b <- function(ny, delta = 1e-09) ifelse(abs(ny) < 
        delta, log(log(2)), log((exp(log(2) * ny) - 1)/ny))
    logz50b <- ny2logz50b(ny, delta)
    ifelse(abs(ny) < delta, log(2)/2, exp(logz50b)/exp(log(1 + 
        ny * exp(logz50b)) * (1/ny + 1)))
}
