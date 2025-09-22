#' Mask bins in CNVisionObj
#'
#' This method masks specified bins in a CNVisionObj object based on given mask types.
#' Masking is used to exclude problematic genomic regions (e.g., centromeres, telomeres)
#' from downstream analyses.
#'
#' @param object A CNVisionObj object containing bin and data information.
#' @param masked_bins An optional integer vector of bin indices that are already masked.
#' @param mask_types A character vector specifying types of regions to mask.
#'                   Allowed values include "none", "centromere", "clone", "contig",
#'                   "heterochromatin", "scaffold", "short_arm", and "telomere".
#'                   Use "none" to disable masking, or NULL to apply default mask types.
#'
#' @return The modified CNVisionObj object with specified bins masked.
#'
#' @examples
#' \dontrun{
#' obj <- Maskbins(obj, mask_types = c("centromere", "telomere"))
#' obj <- Maskbins(obj, mask_types = "none")
#' }
#'
#' @export
setGeneric("Maskbins", function(object,
                                mask_types = c("centromere", "clone",
                                               "contig", "heterochromatin",
                                               "scaffold", "short_arm", "telomere")) {
  standardGeneric("Maskbins")
})

#' @rdname Maskbins
#' @export
setMethod("Maskbins", "CNVisionObj", function(object,
                                              mask_types = c("centromere", "clone",
                                                             "contig", "heterochromatin",
                                                             "scaffold", "short_arm", "telomere")) {
  valid_types <- c("centromere", "clone", "contig", "heterochromatin",
                   "scaffold", "short_arm", "telomere")
  # 1. 判断是否非空（必须至少选 1 个）
  if (length(mask_types) == 0) {
    stop("`mask_types` must include at least one value.")
  }

  # 2. 判断是否合法（是否都属于 valid_types）
  if (!all(mask_types %in% valid_types)) {
    invalid <- mask_types[!mask_types %in% valid_types]
    stop(sprintf("Invalid mask_types: %s\nAllowed values are: %s",
                 paste(invalid, collapse = ", "),
                 paste(valid_types, collapse = ", ")))
  }
  data(gap, package = "CNVision")  # 确保gap数据在包内
  gap <- gap[gap$type %in% mask_types]
  gap_bins <- which(IRanges::countOverlaps(object@bin, gap) > 0)

  object <- FilterBins(object,gap_bins,keep = FALSE)
  return(object)
})


#' Filter bins or rows in CNVisionObj
#' @param object A CNVisionObj object
#' @param index A logical, integer, or character vector indicating bins to keep or drop
#' @param keep If TRUE (default), keep the bins in index; if FALSE, remove them
#' @return A filtered CNVisionObj object
#' @export
FilterBins <- function(object,index,keep = True){
  if(keep){
    object@config$bin_number <- length(index)
  }else{
    object@config$bin_number <- object@config$bin_number -length(index)
    index = -index
  }
  object@cnvData <- object@cnvData[index,, drop = FALSE]
  object@bin <- object@bin[index]
  return(object)
}

#' Filter CNVisionObj by chromosome, position, or cell
#'
#' @param object A CNVisionObj object
#' @param chr A vector of chromosome names (e.g. "chr1", "chrX")
#' @param start,end Optional position range (numeric), used together with `chr`
#' @param cell A vector of sample or cell names to keep/drop (column names in cnvData)
#' @param mode "keep" (default) to retain, or "drop" to exclude matching entries
#' @return A filtered CNVisionObj object
#' @export
FilterBy <- function(object, chr = NULL, start = NULL, end = NULL,
                     cell = NULL, mode = c("keep", "drop")) {
  mode <- match.arg(mode)

  sel_bin <- rep(TRUE, length(object@bin))  # 默认全保留

  # 筛选 chr
  if (!is.null(chr)) {
    sel_chr <- as.character(GenomeInfoDb::seqnames(object@bin)) %in% chr
    sel_bin <- sel_bin & sel_chr
  }

  # 筛选 start/end
  if (!is.null(start) || !is.null(end)) {
    bin_start <- BiocGenerics::start(object@bin)
    bin_end   <- BiocGenerics::end(object@bin)
    if (!is.null(start)) sel_start <- bin_end >= start else sel_start <- TRUE
    if (!is.null(end))   sel_end   <- bin_start <= end  else sel_end <- TRUE
    sel_pos <- sel_start & sel_end
    sel_bin <- sel_bin & sel_pos
  }

  # mode 翻转
  if (mode == "drop") sel_bin <- !sel_bin

  # 调用底层 FilterBins（保留或删除 bin 和对应行）
  #object <- FilterBins(object, index = sel_bin, mode = "keep")

  # 筛选 cell（列过滤）
  if (!is.null(cell)) {
    if (!all(cell %in% object@config$cells)) {
      stop("Some cells not found in cnvData")
    }
    if (mode == "keep") {
      object@cnvData <- object@cnvData[, cell, drop = FALSE]
    } else {
      object@cnvData <- object@cnvData[, !( cell %in% object@config$cells), drop = FALSE]
    }
  }

  return(object)
}
