#' @title Internal: Scale a vector
#' @description Multiply and shift a numeric vector: `x * multiplier + offset`.
#' @param x Numeric vector.
#' @param multiplier Scale factor.
#' @param offset Value to add. 、
#' @return Transformed numeric vector.
#' @noRd
scale_transform <- function(x,object, cell) {
  multiplier <- object@result$optimal_scale[[cell]][,multiplier]
  offset <- object@result$optimal_scale[[cell]][,offset]
  x[,which(object@config$cells == cell)] * multiplier + offset
}
