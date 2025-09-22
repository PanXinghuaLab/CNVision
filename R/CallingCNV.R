#' Calling the CNV
#'
#' @param object  a CNVcaller object
#' @param modelNames A vector of character strings indicating the models to be fitted in the EM phase of clustering.
#'
#' @return object
#' @import mclust
#' @export
#'
#' @examples CallingCNV(object, "E")
CallingCNV <- function(object, modelNames = "E") {
  data <- merge_seg.mean.LOWESS(object)
  # 运行 GMM 聚类
  model <- Mclust(data,
    modelNames = modelNames,
    control = emControl(itmax = 500, tol = 1e-5)
  )
  classification <- model[["classification"]]
  gaussian_params <- data.frame(
    # CN = unique(classification),
    CN = round(tapply(data, classification, mean), 0),
    mean = tapply(data, classification, mean),
    sd = tapply(data, classification, sd),
    weight = as.numeric(table(classification) / length(classification))
  )
  object@options$classification <- gaussian_params
  object@data@CN <- round(object@data@seg.mean.LOWESS * object@options$best_sos_point["multipliers"])
  print(plot_gaussian(object))
  return(object)
}

clip_values <- function(x, lower, upper) {
  x <- pmax(x, lower) # 限制最小值
  x <- pmin(x, upper) # 限制最大值
  return(x)
}

remove_out_of_range <- function(x, lower, upper) {
  x <- x[x >= lower & x <= upper] # 只保留符合条件的值
  return(x)
}

merge_seg.mean.LOWESS <- function(object) {
  data <- as.vector(object@cnvData$seg.mean.LOWESS) * object@config$best_ploidy[["multiplier"]]+object@config$best_ploidy[["offset"]]
  data <- remove_out_of_range(data, lower = 0.5, upper = 5.5)
  return(data)
}

sigma3 <- function() {
  lower_bound <- mu - 3 * sigma
  upper_bound <- mu + 3 * sigma

  outliers <- x[x < lower_bound | x > upper_bound]
  return(outliers)
}
