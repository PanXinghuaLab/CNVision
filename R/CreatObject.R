#' Create a CNVision object
#' @param length The length of reads
#' @param genome the Reference Genome
#' @param dir Directory containing BAM files
#' @param resolution the number of bins
#'
#' @return A CNVision S4 object containing cells information and config
#' @export
#'
#' @importFrom S4Vectors DataFrame
#' @examples CNVision(150)

# 创建 bin 的位置信息并生成 CNVision 对象
CNVision <- function(dir) {
  files <- list.files(dir, pattern = "\\.bam$", full.names = TRUE)
  cells <- tools::file_path_sans_ext(basename(files))
  object <- new("CNVisionObj",
    bin = GenomicRanges::GRanges(),
    cnvData = S4Vectors::DataFrame(),
    result = list(),# 使用空的 DataFrame 作为数据
    config = list(
      dir = dir,
      cells = cells,
      files = files
    )
  )
  return(object)
}


