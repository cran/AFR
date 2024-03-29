% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bg.R
\name{bg}
\alias{bg}
\title{Breusch-Godfrey test [BG test]}
\usage{
bg(
  model,
  order = 1,
  order.by = NULL,
  type = c("Chisq", "F"),
  data = list(),
  fill = 0
)
}
\arguments{
\item{model}{is a (generalized)linear regression model}

\item{order}{integer. maximal order of serial correlation to be tested.}

\item{order.by}{Either a vector z or a formula with a single explanatory variable like ~ z}

\item{type}{the type of test statistic to be returned}

\item{data}{an optional data frame containing the variables in the model}

\item{fill}{starting values for the lagged residuals in the auxiliary regression. By default 0 but can also be set to NA.}
}
\description{
BG test is used to test for autocorrelation in the errors of a regression model
}
\examples{
model <- lm(real_gdp ~ imp + exp + poil + eurkzt + tonia_rate, data = macroKZ)
bg(model)
}
\references{
Mitchel, D. and Zeileis, A. Published 2021-11-07. lmtest package
}
