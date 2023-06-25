#'#' This is some description of this function.
#'@title Using the 7 machine learning algorithms to construct a signature for ICI response
#' @param training Train model and tune parameters (10 times repeated, 5 folds cross-validation, 10 tuneLenth for each parameter)
#' @param validation Compare model performance and pick the best one as the final model
#' @param method A string value to indicate the study type for calculation. Allowed values contain c('scRNAseq', 'bulk_RNAseq').
#' @param sig Number of processors to use when doing the calculations in parallel. This requires to previously load either the parallel or the snow library. If parallel is loaded and this argument is left with its default value (parallel.sz=0) then it will use all available core processors unless we set this argument with a smaller number. If snow is loaded then we must set this argument to a positive integer number that specifies the number of processors to employ in the parallel calculation.

#'
#' @author
#'
#' @return Continuous numerical variables with RCDscore, named by cells or samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("training")
#' data("validation")
#' data("test")
#' data("RCD.Sig_selected")
#' library(parallel)
#' cl <- makePSOCKcluster(16)
#' registerDoParallel(cl)
#' res <- CompareModel.for.ICI.response(training = training,
#'                     validation = validation,
#'                     method = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
#'                     sig = RCD.Sig_selected) # about 30 min
#' stopCluster(cl)
#'
#' print(res[['auc']]) ## ' ' achieves best performance
#' }



CompareModel.for.ICI.response <- function(training, validation, method,sig){

  training <- training[,colnames(training) %in% c('response', sig)]
  validation  <- validation[,colnames(validation) %in% c('response',sig)]

  #7 models adpoted in this study as followings:
  #'nb': navie bayes
  #'svmRadialWeights': Support Vector Machines with Class Weights
  #'rf': random forest
  #'kknn': k-Nearest Neighbors
  #'adaboost':AdaBoost Classification Trees
  #'LogitBoost':Boosted Logistic Regressions
  #'cancerclass': cancerclass


  #Grid search for parameter tuning
  Grid <- list( nb = expand.grid(fL =  c(0,0.5,1,1.5,2.0), usekernel = TRUE, adjust = c(0.5,0.75,1,1.25,1.5)),
                svmRadialWeights = expand.grid(sigma = c(0.0005 ,0.001 ,0.005 ,0.01 ,0.05),C = c( 1 ,3 ,5 ,10 ,20), Weight = c(0.1 ,0.5 ,1 ,2 ,3 ,5 ,10)),
                rf = expand.grid(mtry = c(2,42,83,124,165,205,246,287,328,369)),
                kknn = expand.grid(kmax = c(5,7,9,11,13), distance = 2 , kernel = 'optimal'),
                adaboost = expand.grid(nIter = c(50,100,150,200,250) ,method= c('Adaboost.M1','Real adaboost')),
                LogitBoost = expand.grid(nIter = c(11,21,31,41,51,61,71,81,91,101) )
  )
  TuneLength =  list( nb = nrow(Grid[['nb']]),
                      svmRadialWeights = nrow(Grid[['svmRadialWeights']]) ,
                      rf = nrow(Grid[['rf']]),
                      kknn =nrow(Grid[['kknn']]) ,
                      adaboost = nrow(Grid[['adaboost']]),
                      LogitBoost =  nrow(Grid[['LogitBoost']])
  )


  ##model training with different algorithms
  ls_model <- lapply(method,function(m){
    if(m == 'cancerclass'){ # cancerclass is not avaliable in caret
      pData <- data.frame(class = training$response, sample = rownames(training),row.names = rownames(training))
      phenoData <- new("AnnotatedDataFrame",data=pData)
      Sig.Exp <- t(training[,-1])
      Sig.Exp.train <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
      predictor <- generics::fit(Sig.Exp.train, method = "welch.test")
      model.tune <- predictor
    } else{ # other algorithms were calculated using R package caret

      f = 5  # f folds resampling
      r = 10 # r repeats
      n = f*r

      # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
      seeds <- vector(mode = "list", length = n + 1)
      #the number of tuning parameter
      for(i in 1:n) seeds[[i]] <- sample.int(n=1000, TuneLength[[m]])

      #for the last model
      seeds[[n+1]]<-sample.int(1000, 1)


      ctrl <- trainControl(method="repeatedcv",
                           number = f, ## 5-folds cv
                           summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                           classProbs=TRUE,
                           repeats = r, ## 10-repeats cv,
                           seeds = seeds
      )



      model.tune <- train(response ~ .,
                          data = training,
                          method = m,
                          metric="ROC",
                          trControl=ctrl,
                          tuneGrid = Grid[[m]]
      )
    }
    print(m)
    return(model.tune)
  }
  )

  ##model validation
  auc <- lapply(ls_model,function(model.tune){
    if(class(model.tune) == 'predictor'){
      pData <- data.frame(class = validation$response, sample = rownames(validation),row.names = rownames(validation))
      phenoData <- new("AnnotatedDataFrame",data=pData)
      Sig.Exp <- t(validation[,-1])
      Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
      prediction <- predict(model.tune, Sig.Exp.test,"NR", ngenes=nrow(Sig.Exp), dist = "cor")
      roc <- roc(response  = prediction@prediction[,'class_membership'],
                 predictor = as.numeric(prediction@prediction[,'z'])
      )
      roc_result <- coords(roc, "best")
      auc <- data.frame(ROC=roc$auc, Sens = roc_result$sensitivity, Spec = roc_result$specificity)

    }else {
      prob <- predict(model.tune,validation[,-1],type = "prob")
      pre <- predict(model.tune,validation[,-1])
      test_set <- data.frame(obs = validation$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
      auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
    }

    return(auc)
  }) %>% do.call(rbind,.)

  rownames(auc) <- method

  res <- list()

  res[['model']] <- ls_model
  res[['auc']] <- auc


  return(res)

}


