#' @title Breusch-Pagan test
#' @description
#' Breusch-Pagan test is used to test against heteroskedasticity of a time-series
#' @param model is a (generalized)linear regression model
#' @param varformula a formula describing only the potential explanatory variables for the variance (no dependent variable needed). By default the same explanatory variables are taken as in the main regression model.
#' @param studentize logical. If set to TRUE Koenker's studentized version of the test statistic will be used.
#' @param data an optional data frame containing the variables in the model
#' @import stats
#' @references Torsten, H., Zeileis, A., Farebrother, Richard W., Cummins, C., Millo, G., Mitchell, D., lmtest package
#' Wang, B., 2014, bstats package
#' @export

bp<-function (model, varformula = NULL, studentize = TRUE, data = list())
  {
    dname <- paste(deparse(substitute(model)))
    if (!inherits(model, "model")) {
      X <- if (is.matrix(model$x))
        model$x
      else model.matrix(terms(model), model.frame(model))
      y <- if (is.vector(model$y))
        model$y
      else model.response(model.frame(model))
      Z <- if (is.null(varformula))
        X
      else model.matrix(varformula, data = data)
    }
    else {
      mf <- model.frame(model, data = data)
      y <- model.response(mf)
      X <- model.matrix(model, data = data)
      Z <- if (is.null(varformula))
        X
      else model.matrix(varformula, data = data)
    }
    if (!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in%
                row.names(X))))) {
      allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
      X <- X[allnames, ]
      Z <- Z[allnames, ]
      y <- y[allnames]
    }
    if (ncol(Z) < 2)
      stop("the auxiliary variance regression requires at least an intercept and a regressor")
    k <- ncol(X)
    n <- nrow(X)
    resi <- lm.fit(X, y)$residuals
    sigma2 <- sum(resi^2)/n
    if (studentize) {
      w <- resi^2 - sigma2
      aux <- lm.fit(Z, w)
      bp <- n * sum(aux$fitted.values^2)/sum(w^2)
      method <- "studentized Breusch-Pagan test"
    }
    else {
      f <- resi^2/sigma2 - 1
      aux <- lm.fit(Z, f)
      bp <- 0.5 * sum(aux$fitted.values^2)
      method <- "Breusch-Pagan test"
    }
    message1<-"Homoskedasticity presents. Please use other tests additionally.In case of opposite results study the case further."
    message2<-"Heteroskedasticity presents. Please use others tests additionally.In case of opposite results study the case further."
    names(bp) <- "BP"
    df <- c(df = aux$rank - 1)
    PVAL<-pchisq(bp, df, lower.tail = FALSE)

    if (PVAL>=0.05)
      alternative<-message1
    else
      alternative<-message2
    cat(alternative)

    RVAL <- list(statistic = bp, method = method,
                 p.value = PVAL)
    class(RVAL) <- "htest"
    return(RVAL)

    if (PVAL>=0.05)
      alternative<-message1
    else
      alternative<-message2
    #cat(alternative)

}

