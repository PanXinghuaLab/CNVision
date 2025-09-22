
#' The CNVisionObj Class
#'
#' Define the 'CNVision' class, containing the 'data' class and other slots
#' @slot cnvData DFrame.
#' @slot bin GRanges.
#' @slot config list.
#' @importClassesFrom S4Vectors DFrame
#' @importClassesFrom GenomicRanges GRanges
#' @exportClass CNVisionObj
setClass("CNVisionObj",
         slots = list(
           cnvData = "DFrame",  # Slot for data
           bin = "GRanges",  # Slot for binning information
           result = "list",#  Slot for CNVision result
           config = "list"  # Slot for CNVision config
         )
)

#' Load bins into CNVision object
#'
#' @param object A \code{CNVisionObj} object.
#' @inheritParams LoadGenomeBins
#'
#' @return A \code{CNVisionObj} with updated bin and config slot.
#' @export
#' @rdname LoadBins
setGeneric("LoadBins", function(object, resolution, length, genome) {
  standardGeneric("LoadBins")
})



#' @title Infer ploidy state of the CNVision object
#' @description Estimate sample-wide ploidy using grid search on segment means
#' @param object A CNVisionObj
#' @param ... Optional parameters: m_start, m_end, m_step, b_start, b_end, b_step
#' @return The CNVisionObj with ploidy info stored in `@options`
#' @export
setGeneric("InferPloidy", function(object, ...) {
  standardGeneric("InferPloidy")
})

#' plotCNV generic
#' @param object CNVisionObj
#' @param ... other args
#' @export
setGeneric("plotCNV", function(object, cell,chr = NULL, annotate = TRUE, without_x = TRUE) standardGeneric("plotCNV"))

#' Plot Ploidy Estimation Grid
#'
#' Plot the heatmap of (multiplier, offset) vs error from ploidy estimation.
#' Requires that \code{InferPloidy()} has been run first.
#'
#' @param object A CNVisionObj
#' @return A ggplot2 object (and it will be printed)
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_viridis_c labs theme_minimal
setGeneric("PlotPloidy", function(object,cell) {
  standardGeneric("PlotPloidy")
})


#' Title
#'
#' @param object A CNVisionObj
#' @param cell peaks
#' @param peakHeightThreshold peakHeightThreshold
#' @param kernel kernel
#' @param adjust adjust
#' @param plot plot
#'
#' @returns peaks
#'
#' @examples peaks
setGeneric("detect_peaks", function(object,peakHeightThreshold = 0.05,
                                    minpeakdistance = 0.5,
                                    kernel = "epanechnikov",
                                    adjust = 1,
                                    plot = TRUE) {
  standardGeneric("detect_peaks")
})


#' Title
#'
#' @param object A CNVisionObj
#'
#' @returns model
#'
#' @examples model
setGeneric("laplaceMM", function(object,core_prob = 0.95,dist_type = "gaussian",background_type = "uniform") {
  standardGeneric("laplaceMM")
})


#' Plot CNV Peaks
#'
#' @description
#' This S4 method plots CNV peaks for a \code{CNVisionObj} object.
#' It visualizes both histogram of CNV segments and density/peak information.
#'
#' @param object A \code{CNVisionObj} instance containing CNV analysis results.
#' @param ... Additional arguments passed to plotting functions (optional).
#'
#' @return A \code{ggplot} or \code{patchwork} object displaying CNV peaks.
#'
#' @examples
#' \dontrun{
#'   data(myCNVObj)
#'   plotPeaks(myCNVObj)
#' }
#'
#' @export
#'
setGeneric("plotPeaks", function(object) {
  standardGeneric("plotPeaks")
})
