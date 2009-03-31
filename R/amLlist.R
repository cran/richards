
setClass("amLlist", contains = "list")

`llamodified` <-
function (list, FUN, 
          applyFUN2X = function(X, FUN, ...) FUN(X, ...), ...) 
{
    fits <- lapply(list, function(i)
                lamodified(i, FUN, applyFUN2X, ..., noAttributes = TRUE)
              # lapply(i, function(X) applyFUN2X(X, FUN, ...))
                )
    result <- new("amLlist", fits)
    attr(result, "FUN") <- FUN
    attr(result, "applyFUN2X") <- applyFUN2X
  # attr(result, "Datalist") <- Datalist
    attr(result, "rest") <- list(...)
    return(result)
}

setMethod("coef", "amLlist", function(object, ...) lapply(object, coef) )

setMethod("summary", "amLlist", function(object) lapply(object, summary) )

setMethod("show", "amLlist",
    function(object) { 
        coef <- coef(object)
        cat("Call:\n")
        print(attr(object, "FUN"))
        cat("Coefficients:\n")
        print(coef)
 } )
