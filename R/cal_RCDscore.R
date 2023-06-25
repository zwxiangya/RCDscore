#' #' This is some description of this function.
#' @title  Calculate regulated cell death score using scRNA or bulk transcriptomic data
#' @param expr log2 transformed TPM matrix derived from scRNAseq, or normalized gene expression matrix from bulk samples; cols are genes, rows are cells or samples.
#' @param rcds The list containing selected regulated cell death features.
#' @param study.type A string value to indicate the study type for calculation. Allowed values contain c('scRNAseq', 'bulk_RNAseq').
#' @param parallel.sz Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.
#' @param verbose Gives information about each calculation step. Default: FALSE.

#'
#' @author Wei Zhang
#'
#' @return Continuous numerical variables with RCDscore, named by cells or samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("exprset")
#' data("rcd.sig")
#' RCDscore <- cal_RCDscore(expr = t(exprset), rcds = rcd.sig, study.type = "bulk_RNAseq",parallel.sz = 12,verbose = T)
#' }
cal_RCDscore <- function(expr, # gene expression matrix with genes as colnames and samples or cells as rownames
                         rcds = NULL,
                         study.type = NULL, # "bulk_RNAseq" or "scRNAseq"
                         parallel.sz = NULL,
                         verbose = FALSE) {
  library(GSVA)
  options(warn = -1)
  if (is.null(study.type)) {
    stop("please set your study types as 'scRNAseq' or 'bulk_RNAseq'")
  }
  # prepare input data
  inputMatrix <- expr
  rownames(inputMatrix) <- gsub("-", ".", rownames(inputMatrix))
  bg_genes <- rcds
  # set study types
  # scRNAseq
  if (study.type == "scRNAseq") {
    ssgsea <- gsva(as.matrix(t(inputMatrix)), bg_genes, method = "gsva", kcdf = "Gaussian", parallel.sz = parallel.sz)
  }

  # bulk_RNAseq
  if (study.type == "bulk_RNAseq") {
    ssgsea <- gsva(as.matrix(t(inputMatrix)), bg_genes, method = "ssgsea", kcdf = "Gaussian", parallel.sz = parallel.sz)
  }

  # calculate RCDscore
  ssgsea <- ssgsea %>% as.data.frame()
  RCD.score <- apply(ssgsea, 2, sum) %>% as.data.frame()
  colnames(RCD.score) <- "RCD.score"

  return(RCD.score)
}
