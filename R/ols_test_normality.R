#' Test for normality
#' Test for detecting violation of normality assumption.
#'
#' @param model an object of class \code{lm}.
#' @param ... Other arguments.
#'
#' @return \code{ols_test_normality} is a list containing the
#' following components:
#'
#' \item{kolmogorv}{kolmogorov smirnov statistic}
#' \item{shapiro}{shapiro wilk statistic}
#' \item{cramer}{cramer von mises statistic}
#' \item{anderson}{anderson darling statistic}
#'
#' @examples
#' model <- lm(real_gdp ~ imp + exp + usdkzt + poil, data = macroKZ)
#' ols_test_normality(model)
#'
#' @import olsrr
#' @importFrom stats ks.test shapiro.test
#' @importFrom goftest cvm.test
#' @importFrom nortest ad.test
#'
#' @export
#'
ols_test_normality <- function(model, ...) UseMethod("ols_test_normality")

#' @export
#'
ols_test_normality.default<- function(model, ...) {

  ks  <- ks.test(model, "pnorm", mean(model), sd(model))
  sw  <- shapiro.test(model)
  cvm <- cvm.test(model)
  ad  <- ad.test(model)

  result <- list(kolmogorv = ks,
                 shapiro   = sw,
                 cramer    = cvm,
                 anderson  = ad)

  class(result) <- "ols_test_normality"
  return(result)
}
#' @export
#' @rdname ols_test_normality
#'
ols_test_normality.lm <- function(model, ...) {

  if (!all(class(model) == "lm")) {
    stop("Please specify a OLS linear regression model.", call. = FALSE)
  }

  ols_test_normality.default(residuals(model))

}
#' @export
#'
print.ols_test_normality <- function(x, ...) {
  print_norm_test(x)
}
