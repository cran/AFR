% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gq.R
\name{gq}
\alias{gq}
\title{Godfrey-Quandt test}
\usage{
gq(
  model,
  point = 0.5,
  fraction = 0,
  alternative = c("greater", "two.sided", "less"),
  order.by = NULL,
  data = list()
)
}
\arguments{
\item{model}{is a (generalized)linear regression model}

\item{point}{numerical. If point is smaller than 1 it is interpreted as percentages of data}

\item{fraction}{numerical. The number of central observations to be omitted.}

\item{alternative}{a character string specifying the alternative hypothesis.}

\item{order.by}{Either a vector z or a formula with a single explanatory variable like ~ z}

\item{data}{an optional data frame containing the variables in the model.}
}
\value{
Textual result of the test.
List of values:
\item{GQ}{test statistics}
\item{p-value}{p-value of the BP test}
}
\description{
Godfrey-Quandt test is used to test against heteroskedasticity of a time-series
}
\references{
Torsten, H., Zeileis, A., Farebrother, Richard W., Cummins, C., Millo, G., Mitchell, D., lmtest package
Wang, B., 2014, bstats package
}
