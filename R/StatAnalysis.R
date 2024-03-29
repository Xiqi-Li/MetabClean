
# summary stats

#' Get feature fill rate
#'
#' get percentage of non NA cells per row or per column
#' @param x The data matrix or data frame
#' @param byColumn Boolean: whether the function done by column
#' @return The fill rate
#' @examples
#' fill <- getFill(data_mx, byColumn=F);
#' @export
getFill=function(x,byColumn=F){
  if(byColumn){
    fill=apply(x, 2, function(x) sum(is.na(x)))/nrow(x)
    names(fill)=colnames(x)
  }else{
    fill=apply(x, 1, function(x) sum(is.na(x)))/ncol(x)
    names(fill)=rownames(x)
  }

  return(1-fill)
}

#' Get rate of abnormal measures by column
#'
#' get percentage of abnormal measures per column
#' @param x The data matrix or data frame
#' @param threshold threshold where measures be count as abnormal if its absolute value is above the or below the threshold
#' @return The abnormal rate
#' @examples
#' fill <- getPertubedByColumn(data_mx, threshold=2);
#' @export
getPertubedByColumn=function(x,threshold){
  Nabnormal=apply(x, 2, function(x) sum(na.omit(x)<(-threshold)|na.omit(x)>threshold)/length(na.omit(x)))
  return(Nabnormal)
}

#' Order and cut
#'
#' Decreasingly order absolute values of measures and output top perturbed N features
#' @param x The data matrix or data frame
#' @param threshold threshold where measures be count as abnormal if its absolute value is above the or below the threshold
#' @param N maximum number of features to output
#' @return The abnormal rate
#' @examples
#' top30 <- getPertubedByColumn(data_mx, threshold=2,N=30);
#' @export
orderNCut=function(x,threshold,N){
  x=na.omit(x)
  x=x[abs(x)>threshold]
  x=x[order(abs(x),decreasing = T)]
  if(length(x)>N){x=x[1:N]}
  return(x)
}

#' get association
#'
#' perform profile-wise association analysis
#' @param D The data frame containing predictors, responses and covariates
#' @param ress vector containing all names of all response
#' @param preds vector containing all names of all predictors
#' @param ctrlVs varaibles that need to be controled for
#' @param padjMethod p-value adjustment method
#' @return A data frame of coeffecients and p-values
#' @examples
#' assc = getAssociatio(D, ress=c("age","gender"), preds=grep("feature",colnames(D),value=T));
#' @export
getAssociation=function(D,ress,preds,ctrlVs=NULL,padjMethod="bonferroni"){
  assc=data.frame(outcome=character(),exposure=character(),coef=numeric(),pVal=numeric(),padj=numeric())
  r=0
  for (res in ress){
    for (pred in preds){
      r=r+1
      formu = paste0(res, " ~ ", paste(c(sprintf("`%s`",pred),ctrlVs),collapse = " + "))
      mod=lm(data = D, as.formula(formu))
      assc[r,c("outcome","exposure","coef","pVal")]=c(res,pred,summary(mod)$coefficients[2,1],summary(mod)$coefficients[2,4])
    }
    assc$coef=as.numeric(assc$coef)
    assc$pVal=as.numeric(assc$pVal)
    assc[,"padj"]=p.adjust(assc[,"pVal"],method = padjMethod)#
  }
  return(assc)
}


#' get t-test result
#'
#' perform profile-wise t-test by row
#' @param dataGroup1 The group 1 data matrix. Feature by row, sample by column.
#' @param dataGroup2 The group 2 data matrix. Feature by row, sample by column.
#' @return A data frame of mean differences and p-values
#' @examples
#' t.result = performTtestsAllClassesOneVsRest(dataMatrix=data_mx,classVector=c("treatment1","treatment2","treatment1","control"));
#' @export
performTtestsAllRows = function(dataGroup1,dataGroup2){
  nGroup1 = ncol(dataGroup1)
  nGroup2 = ncol(dataGroup2)
  dataAll = cbind(dataGroup1,dataGroup2)
  tTestWithErrorHandling = function(x){
    if(nGroup1==1){
      testResult = try(t.test(mu=unlist(x[1:nGroup1]),x[(nGroup1+1):(nGroup1+nGroup2)]),silent=TRUE);
      if(is.character(testResult)){
        warning(testResult)
        c(NA,NA,NA)
        }else{
          c(testResult$p.value,unlist(x[1:nGroup1]),testResult$estimate)
        }

      }else if(nGroup2==1){
        testResult = try(t.test(x[1:nGroup1],mu=unlist(x[(nGroup1+1):(nGroup1+nGroup2)])),silent=TRUE);
        if(is.character(testResult)){
          warning(testResult)
          c(NA,NA,NA)
          }else{
            c(testResult$p.value,testResult$estimate,unlist(x[(nGroup1+1):(nGroup1+nGroup2)]))
          }
        }else{
          testResult = try(t.test(x[1:nGroup1],x[(nGroup1+1):(nGroup1+nGroup2)]),silent=TRUE);
          if(is.character(testResult)){
            warning(testResult)
            c(NA,NA,NA)
            }else{
              c(testResult$p.value,testResult$estimate)
            }
        }
  }

  results = matrix(unlist(apply(dataAll,1,tTestWithErrorHandling)),ncol=3,byrow=TRUE)
  colnames(results) = c("P.value","Mean.group.1","Mean.group.2")
  rownames(results) = rownames(dataGroup1)
  results
}

#' get t-test result
#'
#' get t-test result by group, one versus the rest
#' @param dataMatrix The data matrix
#' @param classVector vector matching column names to categories
#' @examples
#' t.result = performTtestsAllClassesOneVsRest(dataMatrix,classVector)
#' @return A data frame of mean differences and p-values
#' @export
performTtestsAllClassesOneVsRest = function(dataMatrix,classVector){
  if(ncol(dataMatrix)!=length(classVector)){
    stop("Number of columns of data matrix must be equal to the length of the class vector")
  }
  possibleClasses = unique(classVector)
  nClasses = length(possibleClasses)

  allPvalues = matrix(NA,nrow=nrow(dataMatrix),ncol=nClasses)
  allDiffMeans = matrix(NA,nrow=nrow(dataMatrix),ncol=nClasses)
  colnames(allPvalues) = possibleClasses
  rownames(allPvalues) = rownames(dataMatrix)
  colnames(allDiffMeans) = possibleClasses
  rownames(allDiffMeans) = rownames(dataMatrix)

  for(i in 1:nClasses){
    class = possibleClasses[i]
    resultTest = performTtestsAllRows(dataMatrix[,classVector==class],dataMatrix[,classVector!=class])
    allPvalues[,i] = resultTest[,1]
    allDiffMeans[,i] = resultTest[,2]-resultTest[,3]
  }
  result = list(allPvalues,allDiffMeans)
  names(result) = c("P.Values","Difference.Between.Means")
  return(result)
}


#' get t-test result by group, pairwise
#'
#' get t-test result by group, pairwise
#' @param dataMatrix The data matrix
#' @param classVector vector matching column names to categories
#' @examples
#' t.result = performTtestsAllClassesEachPair(dataMatrix,classVector)
#' @return A data frame of mean differences and p-values
#' @export
performTtestsAllClassesEachPair = function(dataMatrix,classVector){
  if(ncol(dataMatrix)!=length(classVector)){
    stop("Number of columns of data matrix must be equal to the length of the class vector")
  }
  possibleClasses = unique(classVector)
  nClasses = length(possibleClasses)

  allPValues = NULL
  allDiffMeans = NULL
  names = NULL
  for(i in 1:(nClasses-1)){
    for(j in (i+1):nClasses){
      class1 = possibleClasses[i]
      class2 = possibleClasses[j]
      names = c(names,paste(class1,class2,sep="."))
      result = performTtestsAllRows(dataMatrix[,classVector==class1],dataMatrix[,classVector==class2])
      allPValues = cbind(allPValues,result[,1])
      allDiffMeans = cbind(allDiffMeans,result[,2] - result[,3])
    }
  }
  colnames(allPValues) = names
  colnames(allDiffMeans) = names
  result = list(allPValues,allDiffMeans)
  names(result) = c("P.Values","Difference.Between.Means")
  return(result)
}
