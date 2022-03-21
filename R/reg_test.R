#' Test for detecting violation of Gauss-Markov assumptions.
#'
#' @param y A numeric vector or an object of class \code{lm}.
#'
#' @return \code{reg_test} returns an object of class \code{"reg_test"}.
#' An object of class \code{"reg_test"} is a list containing the
#' following components:
#'
#' \item{bp}{Breusch-Pagan statistic}
#' \item{bg}{Breusch-Godfrey statistic}
#' \item{dw}{Durbin-Watson statistic}
#' \item{gq}{Godfrey-Quandt statistic}

#' @examples
#' model <- lm(real_gdp~imp+exp+poil+eurkzt+usdkzt, macroKZ)
#' reg_test(model)

#' @importFrom lmtest bgtest dwtest
#' @importFrom lmtest bptest gqtest
#'
#' @export
#'


reg_test <- function(y) {


  bp  <- bptest(y)
  bg  <- bgtest(y)
  dw  <- dwtest(y)
  gq  <- gqtest(y)

  result <- list("1"  = bp,
                 "2"  = bg,
                 "3"  = dw,
                 "4"  = gq
                 )

  return(result)
}



