#' @title Data check for errors
#' @description
#' Preliminary check of data frame for missing values, wrong format, outliers.
#' @usage checkdata(x)
#' @param x is a data frame
#' @import stats
#' @importFrom tseries na.remove
#' @importFrom utils data
#' @export

#must be without N/A to run

globalVariables(c("macroKZ"))

checkdata<-
  function (x)

  {
    miss<-sapply(macroKZ, function(x){
      sum(is.na(x))})

    missing_total=0
    for (i in 1:length(miss)){
      missing_total = missing_total+miss[[i]]
    }

    number<-sapply(x, function(x){
      as.logical(sum(is.numeric(x)))})

    out <- function(x){
      count=0
      sdt<-3*sd(x, na.rm=TRUE)
      m<-mean(x)
      for (x in colnames(macroKZ))
      {
        if ((x>m+sdt) & (x<m-sdt))
          count<-count+1
        return(count)}
    }

    outlier<-sapply(na.remove(macroKZ), out)

    outlier_total=0
    for (i in 1:length(outlier)){
      outlier_total = outlier_total+outlier[[i]]
    }

    result <- list("Missing_items" = miss,
                   "Numeric_format"   = number,
                   "Outliers"=outlier)

    notnum_total=0
    for (i in 1:length(number)){
      if (!isTRUE(i))
        number[[i]]=0
      notnum_total = notnum_total+number[[i]]
    }

    cat(paste("There are", missing_total, "missing items in the dataset"),
        paste("There are", notnum_total, "items in non-numeric format in the dataset"),
        paste("There are", outlier_total, "outliers in the dataset"), sep="\n")


    return(result)

  }




