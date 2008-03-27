`initialNyFixedLog` <-
function (mCall, data, LHS, ny = NULL, minnrow = 5, 
          pnames = c("a", "d", "xmid", "scal")) 
{
    xy <- sortedXyData(mCall[["input"]], LHS, data)
    if (!is.numeric(ny)) 
        ny <- eval(mCall$ny)
    if (!is.numeric(ny)) 
        ny <- xy$ny
    if (nrow(xy) < minnrow) {
        stop(paste(
            "too few distinct input values to fit a Richards curve (ny = ", 
            ny, ")", sep = ""))
    }
    fit <- nls(y ~ SSfpl(x, A, B, xmid, scal), 
               control = nls.control(warnOnly = TRUE), data = xy)
    pars <- as.vector(coef(fit))
    .a <- pars[1]
    .d <- pars[2]
    .xmid <- pars[3]
    .b <- solveB(a = .a, d = .d, x50 = exp(.xmid), b4 = 1/pars[4], ny = ny)
    value <- c(.a, .d, .xmid, 1/.b)
    names(value) <- mCall[pnames]
    value
}
