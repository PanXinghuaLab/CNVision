#' Count the reads of each bin
#'
#' @param object a Bin object
#' @param resolution the number of bins
#'
#' @return a CNVcaller object
#' @export
#'
#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom SummarizedExperiment assay
#' @examples
#' CountRead(object)
#' @export
setGeneric("CountRead", function(object, yieldSize = 100000,
                                 mode = "Union",singleEnd = TRUE,
                                 ignore.strand = TRUE,fragments = FALSE,...) {
  standardGeneric("CountRead")
})

#' @rdname CountRead
setMethod("CountRead", "CNVisionObj", function(object,yieldSize = 100000,
                                               mode = "Union",singleEnd = TRUE,
                                               ignore.strand = TRUE,fragments = FALSE,...) {
  bam <- BamFileList(object@config$files, yieldSize = yieldSize,...)
  counts <- summarizeOverlaps(
    features = object@bin,  # 注释区域
    reads = bam,  # BAM 文件
    mode = mode,  # 计数模式，可选 "Union", "IntersectionStrict", "IntersectionNotEmpty"
    singleEnd = singleEnd,  # 单端测序数据
    ignore.strand = ignore.strand,  # 是否忽略链特异性
    fragments = fragments,  # 是否计算片段（用于双端测序）
    ...
  )
  object@cnvData <- DataFrame(row.names = 1:object@config$bin_number)
  object@cnvData$bincount <- assay(counts)

  return(object)
})

