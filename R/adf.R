#' @title Augmented Dickey Fuller Test
#' @description
#' ADF test are used to test stationarity of a time-series data
#' @param x time-series vector
#' @param k the lag order to calculate the test statistic.
#' @import stats
#' @importFrom cli console_width
#' @examples
#' data(macroKZ)
#' adf(macroKZ)
#' @references Trapletti, A., Augmented Dickey-Fuller Test
#' Trapletti, A., KPSS Test for Stationarity
#' @export

adf<-
  function (x,
            k = trunc((length(x) -1)^(1/3)))
  {

    #if ((NCOL(x) > 1) || is.data.frame(x))
      #stop("x is not a vector or univariate time series")
    if (any(is.na(x)))
      stop("NAs in x")
    if (k < 0)
      stop("k negative")
    DNAME <- deparse(substitute(x))
    k <- k + 1
    x <- as.vector(x, mode="double")
    y <- diff(x)
    n <- length(y)
    z <- embed(y, k)
    yt <- z[, 1]
    xt1 <- x[k:n]
    tt <- k:n
    if (k > 1) {
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + 1 + tt + yt1)
    }
    else res <- lm(yt ~ xt1 + 1 + tt)
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
                                                            2]
    table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), c(3.95,
                                                            3.8, 3.73, 3.69, 3.68, 3.66), c(3.6, 3.5, 3.45, 3.43,
                                                                                            3.42, 3.41), c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), c(1.14,
                                                                                                                                                  1.19, 1.22, 1.23, 1.24, 1.25), c(0.8, 0.87, 0.9, 0.92,
                                                                                                                                                                                   0.93, 0.94), c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), c(0.15,
                                                                                                                                                                                                                                        0.24, 0.28, 0.31, 0.32, 0.33))
    table <- -table
    tablen <- dim(table)[2]
    tableT <- c(25, 50, 100, 250, 500, 1e+05)
    tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
    tableipl <- numeric(tablen)
    for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[,
                                                              i], n, rule = 2)$y
    interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
    PVAL<-interpol

    #print
    a <- c("Dickey-Fuller", "p-value")
    b <- c(round(STAT, 3), round(PVAL,3))
    w1 <- max(nchar(a))
    w2 <- max(nchar(b))
    w3 <- console_width()
    w <- sum(w1, w2, 7)
    n <- length(b)


    cat(format(as.character("Augmented Dickey-Fuller test"), width=w3, justify="centre"), "\n")
    if (PVAL>=0.1)
      cat(paste("The result is non-stationary.", "\n"))
    else
      cat(paste("The result is stationary.",  "\n"))

    cat(rep("-", w), sep = "", "\n")
    for (i in seq(n)) {
      cat(fl(a[i], w1),fsp(),fsp(), fg(b[i], w2), "\n")
    }
    cat(rep("-", w), sep = "", "\n")

}

#adf<-function(x){

  #sapply(x, adf_x)

  #s=0
    #for (i in 1:length(adf)){
      #if (alternative==message1)
        #s = s+1
    #}

  #ns=0
    #for (i in 1:length(adf)){
      #if (alternative==message2)
        #ns = ns+1
#}

    #cat(paste("There are", s, "stationary items in the dataset"),
        #paste("There are",ns, "non-stationary items"),sep="\n")

    #structure(list(p.value = PVAL, result = alternative))

    #result <- list(p.value = PVAL,
                   #result   = alternative)

    #return(result)
  #}

