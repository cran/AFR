#' @title Transforming time-series data to stationary
#' @description
#' Percent change is a change between two consecutive terms, %
#' @usage pct1(x)
#' @import stats
#' @param x time-series vector(s)
#' @rdname pct1
#' @export

pct1<-function(x)((x/stats::lag(x)-1)*100)

#' @title Transforming time-series data to stationary
#' @description
#' Percent change is a change between a term and its lagged value for prior period, %
#' @usage pct4(x)
#' @import stats
#' @param x time-series vector(s)
#' @rdname pct4
#' @export

pct4<-function(x){
  x4<-stats::lag(x,4)
  p<-(x/x4-1)*100
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

#' @title Autocorrelation function
#' @description
#' Autocorrelation function demonstrates correlation between stationary time series and its lagged values
#' @param x time-series vector
#' @param lag.max maximum lag at which to calculate the ACF
#' @param type character string giving the type of ACF to be computed
#' @param plot logical
#' @param na.action function to be called to handle missing values. na.pass can be used
#' @param demean logical
#' @param ... other arguments
#' @references Gilbert, P.,  Plummer, M., Ripley, B.D., Auto- and Cross- Covariance and -Correlation Function Estimation
#' Dancho, Matt. Published 2021-01-18. timetk package
#' @import stats
#' @rdname ACF
#' @export

globalVariables(c("C_acf", "plot.acf"))

ACF<-function (x, lag.max = NULL, type = c("correlation", "covariance",
                                           "partial"), plot = TRUE, na.action = na.contiguous, demean = TRUE,
               ...)
{
  cat("Hint: Spikes near zero range demostrate not statistical significance and absence of high correlation")
  type <- match.arg(type)
  if (type == "partial") {
    m <- match.call()
    m[[1L]] <- quote(stats::pacf)
    m$type <- NULL
    return(eval(m, parent.frame()))
  }
  series <- deparse1(substitute(x))
  x <- na.action(as.ts(x))
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("'x' must be numeric")
  sampleT <- as.integer(nrow(x))
  nser <- as.integer(ncol(x))
  if (is.na(sampleT) || is.na(nser))
    stop("'sampleT' and 'nser' must be integer")
  if (is.null(lag.max))
    lag.max <- floor(10 * (log10(sampleT) - log10(nser)))
  lag.max <- as.integer(min(lag.max, sampleT - 1L))
  if (is.na(lag.max) || lag.max < 0)
    stop("'lag.max' must be at least 0")
  if (demean)
    x <- sweep(x, 2, colMeans(x, na.rm = TRUE), check.margin = FALSE)
  lag <- matrix(1, nser, nser)
  lag[lower.tri(lag)] <- -1
  acf <- .Call(C_acf, x, lag.max, type == "correlation")
  lag <- outer(0:lag.max, lag/x.freq)
  acf.out <- structure(list(acf = acf, type = type, n.used = sampleT,
                            lag = lag, series = series, snames = colnames(x)), class = "acf")
  if (plot) {
    plot.acf(acf.out, ...)
    invisible(acf.out)
  }
  else acf.out
}


#' @title Partial autocorrelation function
#' @description
#' (partial) Autocorrelation function demonstrates (partial) correlation between stationary time series and its lagged values
#' @param x time-series vector
#' @param lag.max maximum lag at which to calculate the ACF
#' @param plot logical
#' @param na.action function to be called to handle missing values. na.pass can be used
#' @param demean logical
#' @param ... other arguments
#' @references Gilbert, P.,  Plummer, M., Ripley, B.D., Auto- and Cross- Covariance and -Correlation Function Estimation
#' Dancho, Matt. Published 2021-01-18. timetk package
#' @import stats
#' @rdname PACF
#' @export

globalVariables(c("seasonalaxis"))

PACF<-
  function (x, lag.max = NULL, plot = TRUE, na.action = na.contiguous,
            demean = TRUE, ...)
  {
    cat("Hint: Spikes near zero range demostrate not statistical significance and absence of high correlation")
    object <- acf(x, lag.max = lag.max, type = "partial", na.action = na.action,
                  demean = demean, plot = FALSE)
    object$series <- deparse(substitute(x))
    if (plot) {
      nlags <- dim(object$lag)[1]
      plot.out <- object
      input_list <- as.list(substitute(list(...)))
      ylimarg <- is.element("ylim", names(input_list))
      if (ylimarg) {
        plot(plot.out, xaxt = "n", ...)
      }
      else {
        ylim <- c(-1, 1) * 3/sqrt(length(x))
        ylim <- range(ylim, plot.out$acf)
        plot(plot.out, ylim = ylim, xaxt = "n", ...)
      }
      if (is.element("msts", class(x))) {
        seasonalaxis(attributes(x)$msts, nlags, type = "acf")
      }
      else {
        seasonalaxis(frequency(x), nlags, type = "acf")
      }
      return(invisible(object))
    }
    else {
      return(object)
    }
  }


#' @title Hodrick-Prescott filter
#' @description
#' Hodrick-Prescott filter is a data smoothing technique that removes trending in time series data frame
#' @param x time-series vector
#' @param type character, indicating the filter type
#' @param freq integer
#' @param drift logical
#' @import stats
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

