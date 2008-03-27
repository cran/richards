`initialRichards` <-
function (mCall, data, LHS, ny = 1, minnrow = 6, factor = 1.35, 
    pnames = c("a", "d", "b", "x50", "ny")) 
{
    xy <- sortedXyData(mCall[["input"]], LHS, data)
    if (nrow(xy) < minnrow) {
        stop("too few distinct input values to fit a Richards function")
    }
    xy <- cbind(xy, logx = log(xy$x)/factor)
    fit <- nls(y ~ SSfpl(logx, A, B, xmid, scal), 
               control = nls.control(warnOnly = TRUE), data = xy)
    pars <- as.vector(coef(fit))
    value <- c(pars[1], pars[2], 1/pars[4], exp(pars[3]), ny)
    names(value) <- mCall[pnames]
    value
}
