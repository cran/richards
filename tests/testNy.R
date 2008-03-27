
library(richards)

# source("sourceModified.R")

demo(sourceRichards)

demo(sourceModified)

# source("../demo/sourceNyExtras.R")

demo(sourceNyExtras)


ls()

show <- function(ny, a.1, a.2, label, rlim = 10^-13, dlim = 10^-14) {

  print(label)

  print(c(a.1))
  print(c(a.2))

  print(c(a.1 - a.2))

  if (TRUE) {
    d1 <- attributes(a.1)$gradient - attributes(a.2)$gradient
    s1 <- attributes(a.1)$gradient + attributes(a.2)$gradient
    r1 <- ifelse(d1 == 0, 0, d1 / s1)
    r1 <- ifelse(abs(r1) < rlim, 0, r1)
    if ((any(range(r1, na.rm = TRUE)) != 0) & (any(abs(range(d1, na.rm = TRUE)) > dlim))) {
      print(paste("1/3 -", label, ": ", paste(range(r1, na.rm = TRUE), collapse = ";"),
                                 " + ", paste(range(d1, na.rm = TRUE), collapse = ";")
           ))
    # print(abs(range(d1)))
    # print(r1)
    # print(a.1)
    # print(a.2)
    }
  }

}

fewTests <-  function(a = 0.1, d = 2.4, x50 = 100, b = 2, ny = -0.25,
                      x = c(25, 50, 100, 200, 400, 800)) {

  subtitle <- paste("a =", a, ", d =",   d, ", x50 =", x50, 
                  ", b =", b, ", and ny =",  ny, ".")

  print(subtitle)

  s.2.1 <-        SSnym1(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       SSny.m1(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Demo: SSny.m1")

  s.2.1 <-        SSnym1(     x, a, d, x50 = x50, b = b)
  s.2.2 <-     SS.d.nym1(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Documentation: SS.d.nym1")

  s.2.1 <-        SSnym1(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       sSny.m1(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Simple: sSny.m1")

  s.2.1 <-        SSnym1o3(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       SSny.m1o3(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Demo: SSny.m1o3")

  s.2.1 <-        SSnym1o3(     x, a, d, x50 = x50, b = b)
  s.2.2 <-     SS.d.nym1o3(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Documentation: SS.d.nym1o3")

  s.2.1 <-        SSnym1o3(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       sSny.m1o3(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Simple: sSny.m1o3")


  s.2.1 <-        SSny1(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       SSny.1(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Demo: SSny.1")

  s.2.1 <-        SSny1(     x, a, d, x50 = x50, b = b)
  s.2.2 <-     SS.d.ny1(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Documentation: SS.d.ny1")

  s.2.1 <-        SSny1(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       sSny.1(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Simple: sSny.1")

  s.2.1 <-        SSny2(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       SSny.2(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Demo: SSny.2")

  s.2.1 <-        SSny2(     x, a, d, x50 = x50, b = b)
  s.2.2 <-     SS.d.ny2(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Documentation: SS.d.ny2")

  s.2.1 <-        SSny2(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       sSny.2(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Simple: sSny.2")

  s.2.1 <-        SSny3(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       SSny.3(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Demo: SSny.3")

  s.2.1 <-        SSny3(     x, a, d, x50 = x50, b = b)
  s.2.2 <-     SS.d.ny3(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Documentation: SS.d.ny3")

  s.2.1 <-        SSny3(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       sSny.3(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Simple: sSny.3")

  s.2.1 <-        SSny4(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       SSny.4(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Demo: SSny.4")

  s.2.1 <-        SSny4(     x, a, d, x50 = x50, b = b)
  s.2.2 <-     SS.d.ny4(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Documentation: SS.d.ny4")

  s.2.1 <-        SSny4(     x, a, d, x50 = x50, b = b)
  s.2.2 <-       sSny.4(     x, a, d, x50 = x50, b = b)
  show(ny, s.2.1, s.2.2, "Simple: sSny.4")

}

fewTests()

q()
