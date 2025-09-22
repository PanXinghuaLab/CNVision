#' @rdname NormalizeData
NormalizeByMeanPlus1 <- function(x) {
  x <- x + 1
  x <- x / mean(x)
}

#' Correct Data
#'
#' @param object  a CNVcller object
#' @param method the normalize method
#'
#' @return object
#' @export
#'
#' @examples
#' NormalizeData(object)
NormalizeData <- function(object, method = "Normalize"){
  if (method == "Normalize") {
  object@cnvData$ratio <- apply(object@cnvData$bincount, 2, NormalizeByMeanPlus1)
  }
  gc.ratio <- mapply(lowess.gc,
    x = list(object@bin$gc.content),
    y = split(object@cnvData$ratio, col(object@cnvData$ratio))
  )
  object@cnvData$gc.ratio <- gc.ratio
  return(object)
}

lowess.gc <- function(x, y) {
  low <- lowess(x, log(y), f = 0.05)
  z <- approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}
