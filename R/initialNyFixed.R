`initialNyFixed` <-
function (mCall, data, LHS, ny = NULL, minnrow = 5, 
          pnames = c("a", "d", "b", "x50")) 
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
    xy <- cbind(xy, logx = log(xy$x))
    fit <- nls(y ~ SSfpl(logx, A, B, xmid, scal), 
               control = nls.control(warnOnly = TRUE), data = xy)
    pars <- as.vector(coef(fit))
    .a <- pars[1]
    .d <- pars[2]
    .x50 <- exp(pars[3])
    .b <- solveB(a = .a, d = .d, x50 = .x50, b4 = 1/pars[4], ny = ny)
    value <- c(.a, .d, .b, .x50)
    names(value) <- mCall[pnames]
    value
}
