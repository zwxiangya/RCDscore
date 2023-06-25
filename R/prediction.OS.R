#'#' This is some description of this function.
#'@title Calculation of the risk score of the overall survivla
#' @param exp.matrix A cohort matrix with the genes as the colnames, sample names as the rownamses. The colnames inlcudes the RCD.Sur.Sig, The top 3 colnames are "ID", "OS.time", "OS"
#' @author Wei Zhang
#'
#' @return A dataframe that contains the risk score
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("fit")
#' data("rid")
#' data("exp_test")
#' dd = prediction.OS(exp.matrix = exp_test)
#'
#' }





prediction.OS<- function(exp.matrix=NULL){
  library(dplyr)
  mm = list(exp = exp.matrix)
  mm <- lapply(mm,function(x){
    colnames(x) <- gsub("-",".",colnames(x))
    return(x)})

  mm <- lapply(mm,function(x){
    x[,-c(1:3)] <- apply(x[,-c(1:3)],2,as.numeric)
    return(x)})

  mm <- lapply(mm,function(x){
    x[,c(1:2)] <- apply(x[,c(1:2)],2,as.factor)
    return(x)})

  mm <- lapply(mm,function(x){
    x[,c(2:3)] <- apply(x[,c(2:3)],2,as.numeric)
    return(x)})

  mm <- lapply(mm,function(x){
    x[,-c(1:3)] =  x[,-c(1:3)]  %>%
      mutate_all(~if_else(is.na(.), mean(., na.rm = TRUE), .))
    x = na.omit(x)
    return(x)})

  mm <- lapply(mm,function(x){

    if(length(intersect(rid,colnames(x))) == length(rid)){

      y= x
    } else {
      ls = rid[!rid%in%colnames(x)]
      a = matrix(data = rep(0, length(ls)*nrow(x)), nrow = nrow(x), ncol = length(ls)) %>% as.data.frame()
      colnames(a) = ls
      y = cbind(x,a)
    }


    return(y)}
  )

  mm <- lapply(mm, function(x){x[, c("ID",'OS.time', 'OS', rid)]})
  rs <- lapply(mm, function(x){cbind(x[, 1:3], RS = stats::predict(fit, newdata = x)$predicted)})
  rs.table = rs$exp
  return(rs.table)


}










