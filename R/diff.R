#' @title Transforming time-series data to stationary
#' @description
#' Percent change is a change between two consecutive terms, %
#' @usage pct1(x)
#' @import stats
#' @param x time-series vector(s)
#' @examples
#' data (macroKZ)
#' new<-pct1(macroKZ)
#' @rdname pct1
#' @export

pct1<-function(x)((x/stats::lag(x)-1)*100)

#' @title Transforming time-series data to stationary
#' @description
#' Percent change is a change between a term and its lagged value for prior period, %
#' @usage pct4(x)
#' @import stats
#' @param x time-series vector(s)
#' @examples
#' data (macroKZ)
#' new<-pct4(macroKZ)
#' @rdname pct4
#' @export

pct4<-function(x){
  x4<-stats::lag(x,4)
  p<-(x4/x-1)*100
  return(p)
}

#tr<-function(x){
  #trend<-rollaply(m,width=period,fill=NA, align="center", FUN=mean,na.rm=TRUE)
  #season<-m-trend
  #figure<-numeric(period)
  #l<-length(m)
  #index<-seq.int(1,l,by=period)-1
  #for (i in 1:period) figure[i]<-median(season)
#}

#' @title Transforming time-series data to stationary
#' @description
#' Difference of logarithms is finding the difference between two consecutive logarithm values of a time-series
#' @param x time-series vector
#' @param difference difference between x items
#' @param lag lagged period
#' @import stats
#' @importFrom rlang abort
#' @importFrom xts diff.xts
#' @examples
#' data (macroKZ)
#' new<-pct1(macroKZ)
#' @rdname difflog
#' @export

difflog<-
  function (x, lag = 1, difference = 1)
  {
    if (!is.numeric(x))
      rlang::abort("Non-numeric data detected. 'x' must be numeric.")
    x<-log(x)
    ret_vec <- xts::diff.xts(x = x, lag = lag, differences = difference,
                             arithmetic = TRUE, na.pad = TRUE)
    pad_len <- length(x) - length(ret_vec)
    if (pad_len > 0) {
      ret_vec <- c(rep(NA, pad_len), ret_vec)
    }
    return(ret_vec)
  }

#' @title Hodrick-Prescott filter
#' @description
#' Hodrick-Prescott filter is a data smoothing technique that removes trending in time series data frame
#' @param x time-series vector
#' @param type character, indicating the filter type
#' @param freq integer
#' @param drift logical
#' @import stats
#' @examples
#' data(macroKZ)
#' HP(macroKZ[,2])
#' @rdname HP
#' @export

globalVariables(c("undrift"))

HP<-
  function (x, freq = NULL, type = c("lambda", "frequency"), drift = FALSE)
  {
    if (is.null(drift))
      drift <- FALSE
    xname = deparse(substitute(x))
    type = match.arg(type)
    if (is.null(type))
      type <- "lambda"
    if (is.ts(x)) {
      tsp.x <- tsp(x)
      frq.x <- frequency(x)
      if (type == "lambda") {
        if (is.null(freq)) {
          if (frq.x == 1)
            lambda = 6
          if (frq.x == 4)
            lambda = 1600
          if (frq.x == 12)
            lambda = 129600
        }
        else lambda = freq
      }
    }
    else {
      if (type == "lambda") {
        if (is.null(freq))
          stop("freq is NULL")
        else lambda = freq
      }
    }
    if (type == "frequency") {
      if (is.null(freq))
        stop("freq is NULL")
      else lambda = (2 * sin(pi/freq))^-4
    }
    xo = x
    x = as.matrix(x)
    if (drift)
      x = undrift(x)
    n = length(x)
    imat = diag(n)
    Ln = rbind(matrix(0, 1, n), diag(1, n - 1, n))
    Ln = (imat - Ln) %*% (imat - Ln)
    Q = t(Ln[3:n, ])
    SIGMA.R = t(Q) %*% Q
    SIGMA.n = diag(n - 2)
    g = t(Q) %*% as.matrix(x)
    b = solve(SIGMA.n + lambda * SIGMA.R, g)
    x.cycle = c(lambda * Q %*% b)
    x.trend = x - x.cycle
    if (is.ts(xo)) {
      tsp.x = tsp(xo)
      x.cycle = ts(x.cycle, start = tsp.x[1], frequency = tsp.x[3])
      x.trend = ts(x.trend, start = tsp.x[1], frequency = tsp.x[3])
      x = ts(x, start = tsp.x[1], frequency = tsp.x[3])
    }
    A = lambda * Q %*% solve(SIGMA.n + lambda * SIGMA.R) %*%
      t(Q)
    res <- list(cycle = x.cycle, trend = x.trend, fmatrix = A,
                title = "Hodrick-Prescott Filter", xname = xname, call = as.call(match.call()),
                type = type, lambda = lambda, method = "hpfilter", x = x)
    return(structure(res, class = "mFilter"))
  }

