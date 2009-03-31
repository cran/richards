
# .onLoad <- function(lib, pkg){
#     require(methods)
# }

setClass("amList", contains = "list")

`lamodified` <-
function (list, FUN,
          applyFUN2X = function(X, FUN, ...) FUN(X, ...), ...) 
{
    fits <- lapply(list, function(X) applyFUN2X(X, FUN, ...))
    result <- new("amList", fits)
    if (!is.element("noAttributes", names(list(...)))) {
      attr(result, "FUN") <- FUN
      attr(result, "applyFUN2X") <- applyFUN2X
    # attr(result, "Datalist") <- Datalist
      attr(result, "rest") <- list(...)
    }
    return(result)
}

if (!isGeneric("coef")) {
    setGeneric("coef", function(object, ...)
        standardGeneric("coef")) } 

setMethod("coef", "amList",
    function(object, ...) { 
        invalidFit <- function(fit)
           is.null(fit) | any(class(fit) == c("try-error", "data-error"))
        lengths <- lapply(object, function(fit)
                       if (invalidFit(fit)) 0 else length(coef(fit)))
        N <- max(unlist(lengths))

      # X <- NULL
      # lapply(object, function(fit) {
      #                               if (invalidFit(fit))
      #                                 X <<- rbind(X, rep(NA, N))
      #                               else
      #                                 X <<- rbind(X, z = coef(fit))
      #                               NULL } )
      # dimnames(X)[[1]] <- names(fits)

        names <- lapply(object, function(fit)
                       if (invalidFit(fit)) NULL else names(coef(fit)))
        names <- unique(unlist(names))
        x <- rep(NA, length(names))
        names(x) <- names
        X <- lapply(object, function(fit) {
            if (invalidFit(fit)) x
            else { z <- x; y <- coef(fit); z[names(y)] <- y; z}})
        X <- matrix(unlist(X), ncol = N, byrow = TRUE)
        dimnames(X) <- list(names(object), names)
        return(X)
 } )


if (!isGeneric("summary")) {
    setGeneric("summary", function(object, ...)
        standardGeneric("summary")) } 

setMethod("summary", "amList",
    function(object, ...) { 
        invalidFit <- function(fit)
           is.null(fit) | any(class(fit) == c("try-error", "data-error"))

        lengths <- lapply(object, function(fit)
                       if (invalidFit(fit)) 0 else length(coef(fit)))
        N <- max(unlist(lengths))
        min <- min(which(lengths > 0))

        names <- lapply(object, function(fit)
                       if (invalidFit(fit)) NULL else names(coef(fit)))
        names <- unique(unlist(names))
        x <- rep(NA, length(names))
        names(x) <- names
        X <- lapply(object, function(fit) {
            if (invalidFit(fit)) NULL
            else { summary(fit)$parameters}})
        result <- NULL
        for (name in names) {
            Y <- lapply(X, function(x) 
            if (is.null(x))
                rep(NA, 4)
            else
               x[dimnames(x)[1][[1]] == name,])
            Z <- matrix(unlist(Y), ncol = 4, byrow = TRUE)
            dimnames(Z) <- list(names(object), names(Y[[min]]))
            result <- append(result, list(Z))
        }
        names(result) <- names
        return(list(call = summary(object[[1]])$formula,
                    coefs = result))
} )

if (!isGeneric("show")) {
    setGeneric("show", function(object)
        standardGeneric("show")) } 

setMethod("show", "amList",
    function(object) { 
        invalidFit <- function(fit)
           is.null(fit) | any(class(fit) == c("try-error", "data-error"))

        cat("Call:\n")
        cat("    ")
        lengths <- lapply(object, function(fit)
                       if (invalidFit(fit)) 0 else length(coef(fit)))
        print(summary(object[[min(which(lengths > 0))]])$formula)

        cat("Coefficients:\n")
        coef <- coef(object)
        print(coef)
 } )


