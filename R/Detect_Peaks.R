
setMethod("detect_peaks", "CNVisionObj", function(object,peakHeightThreshold = 0.05,
                                                  minpeakdistance = 0.5,
                                                  kernel = "epanechnikov",
                                                  adjust = 1,
                                                  plot = TRUE) {
  data <- do.call(rbind,object@result$ScaledCNV)
  df <- data[data$seg.mean>= 0.5 & data$seg.mean <= 6.5,]

  result <- .detect_peaks(df,adjust = adjust,kernel = kernel,peakHeightThreshold  = peakHeightThreshold,plot = plot,minpeakdistance= minpeakdistance)
  object@result$kde_peaks <- result$kde_peaks
  object@result$offpeaks <- result$offpeaks
  object@result$dens <- result$dens
  object@result$df <- df

  return(object)
})


#' Title
#'
#' @param df data
#' @param peakHeightThreshold peakHeightThreshold
#' @param kernel kernel
#' @param adjust adjust
#' @param plot TRUE
#'
#' @returns  kde_peaks
#' @export
#'
#' @examples .detect_peaks(x,adjust=1,kernel = "epanechnikov",peakHeightThreshold  = 0.01)
.detect_peaks <- function(df,
                         peakHeightThreshold = 0.05,
                         minpeakdistance ,
                         kernel = "epanechnikov",
                         adjust = 1,
                         plot = TRUE) {
  dens <- density(x = df$seg.mean,weights = df$num.mark/sum(df$num.mark),adjust = adjust, kernel = kernel,n = 1024,bw = 0.1)
  peaks_mat <- pracma::findpeaks(dens$y,zero = "+",
                                 minpeakdistance = 1024/diff(range(df$seg.mean))*minpeakdistance,
                                 minpeakheight = max(dens$y)* peakHeightThreshold)
                                 #minpeakheight = quantile(dens$y, peakHeightThreshold))

  if (!is.null(peaks_mat)) {
    peak_idx <- peaks_mat[, 2]
    left_idx <- peaks_mat[, 3]
    right_idx <- peaks_mat[, 4]

    # 构建数据框
    kde_peaks <- data.frame(
      peak_height = peaks_mat[, 1],
      peak_x = dens$x[peak_idx],
      left_base = dens$x[left_idx],
      right_base = dens$x[right_idx],
      left_y = dens$y[left_idx],
      right_y = dens$y[right_idx]
    )
  }
  kde_peaks <- kde_peaks[order(kde_peaks[[2]], decreasing = FALSE), ]

  pick_nearest_by_round <- function(x, tol = 0.3) {
    ints <- sort(unique(round(x)))
    idx <- sapply(ints, function(t) {
      if (length(x) == 0) return(NA_integer_)
      m <- which.min(abs(x - t))
      if (abs(x[m] - t) <= tol) {
        res <- m
        x[m] <<- Inf   # 占位，避免重复
        return(res)
      } else {
        return(NA_integer_)
      }
    })
    names(idx) <- ints
    idx <- na.omit(idx)
    return(idx)
  }
  idx<-  pick_nearest_by_round(kde_peaks$peak_x)
  offpeaks <- kde_peaks[!(seq_along(kde_peaks$peak_x) %in% idx),]
  kde_peaks <- kde_peaks[idx,]

  if (plot) .PlotPeaks(series = rep(df$seg.mean, times = df$num.mark), dens, kde_peaks,offpeaks)

  return(list(kde_peaks = kde_peaks,offpeaks = offpeaks,dens=dens))
}
