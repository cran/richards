`initialRichardsLogB` <-
function (mCall, data, LHS, ny = 1, minnrow = 6, 
          pnames = c("a", "d", "xmid", "b", "ny")) 
{
    xy <- sortedXyData(mCall[["input"]], LHS, data)
    if (nrow(xy) < minnrow) {
        stop("too few distinct input values to fit a Richards function")
    }
    fit <- nls(y ~ SSfpl(x, A, B, xmid, scal), 
               control = nls.control(warnOnly = TRUE), data = xy)
    pars <- as.vector(coef(fit))
    value <- c(pars[1], pars[2], pars[3], 1/pars[4], ny)
    names(value) <- mCall[pnames]
    value
}
