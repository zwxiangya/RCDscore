#'#' This is some description of this function.
#'@title Testing the ICI model in the testing data and the individual independent datasets
#' @param exp.matrix A cohort matrix with the genes as the colnames, sample names as the rownamses. The colnames inlcudes the RCD.Sig
#' @param parallel.sz Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.
#' @author
#'
#' @return A dataframe that contains the predicted response to the ICI
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("test")
#' data("RCD.Sig_selected")
#' dd = prediction.ICI.response(exp.matrix = test,parallel.sz=16)
#' }





prediction.ICI.response<- function(exp.matrix=NULL,
                                   parallel.sz= NULL
                                   ){
  data("training")
  data("validation")
  data("RCD.Sig_selected")
  library(caret)
  library(parallel)
  library(doParallel)
  cl <- makePSOCKcluster(parallel.sz)
  registerDoParallel(cl)
  res <- CompareModel.for.ICI.response(training = training,
                                       validation = validation,
                                       # method = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                       method = c('svmRadialWeights'),
                                       sig = RCD.Sig_selected)
  stopCluster(cl)

  new=exp.matrix
  prob <- predict(res[['model']],new[,-1],type = "prob")
  pre <- predict(res[['model']],new[,-1])
  test_set <- data.frame(R = prob[[1]][,'R'], pred=pre[[1]])
  rownames(test_set) = rownames(new)
  return(test_set)
}










