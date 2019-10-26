#' Cross-validation for plbpsm
#'
#' Performs k-fold cross validation for GGAMs
#'
#' @importFrom caret createFolds
#' @importFrom stats predict
#' @param formula a PLBPSM formula. These are exactly like the formula for
#' a GLM except that univariable and bivariate smooth terms.
#' @param data A data frame or list containing the model response variable and
#' covariates required by the formula. By default the variables are taken from
#' environment(\code{formula}).
#' @param kfold number of folds splitted
#' @param repeat.times repeated times for cross validation
#' @param group whether model identification is conducted for each covariates
#' @param family This is a family object specifying the distribution and link to use in fitting
#' etc. See glm and family for more details.
#' @param ... other arguments.
#' @export
cv.plbpsm<-function(formula, data, kfold = 10,#type.measure = c("mse","deviance", "class", "auc", "mae"),
                    repeat.times=1,family,group=NULL,...)
  {
  n=dim(data)[1]
  cverror_all=c()
  for (i in 1:repeat.times){
    #set.seed(2*i)
    flds <- createFolds(1:n, k = kfold, list = TRUE, returnTrain = FALSE)
    pred=rep(NA,n)
    for(kk in 1:length(flds)){
      dat.est=data[-flds[[kk]],]
      res_gauss2=plbpsm(formula=formula,data=as.data.frame(dat.est),family=family,group=group)
      # a=predict(res_gauss2,as.data.frame(dat_gauss_GGAMS[c(((f-1)*100+1):(f*10)),]))
      pred[flds[[kk]]]=predict(res_gauss2,as.data.frame(data[flds[[kk]],]))
    }
    cverror=sum((pred-data[,as.character(formula[[2]])] )^2,na.rm = TRUE)/sum(!is.na(pred))
    cverror_all=c(cverror_all,cverror)
    result=list(cverror=mean(cverror_all),pred=pred,foldind=flds)
  }
   return(result)
}
