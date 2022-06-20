
makeLegalRowname = function(mets){
  rown = tolower(trimws(mets))
  rown = sapply(rown, function(x) sub("\\*","",x, perl = FALSE))
  rown = trimws(rown)
  return(rown)
}

cleanUpRownames = function(mat, useRownames=FALSE) {
  if (useRownames) {
    tmp = rownames(mat)
  } else {
    tmp = mat[,1]
  }
  tmp = gsub("\\*", "", tmp)
  tmp = gsub("\\\"", "", tmp)
  tmp = gsub(" \\(1)", "", tmp)
  tmp = gsub(" \\(2)", "", tmp)
  tmp = gsub(" \\(3)", "", tmp)
  tmp = gsub(" \\[1]", "", tmp)
  tmp = gsub(" \\[2]", "", tmp)
  ind = which(duplicated(tmp))
  if (length(ind)>0) {
    mat = as.matrix(mat[-ind,])
    tmp = tmp[-ind]
  }
  rownames(mat) = as.character(tolower(sapply(tmp, trimws)))
  if (!useRownames) {
    mat = as.matrix(mat[,-1])
  }
  return(mat)
}

imputeMissingValues = function(data, ref) {
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = data[sort(rownames(data)),]
  ref = ref[sort(rownames(ref)),]
  imputed.data = data
  for (met in 1:nrow(ref)) {
    rowData = ref[met,]
    if (any(is.na(rowData))) {
      rowData = as.numeric(rowData[-which(is.na(rowData))])
    } else {
      rowData = as.numeric(rowData)
    }
    # Impute using uniform random variable, where a = 0.99*observed minimum, and b = observed minimum
    min_row = min(rowData)
    if (min_row<0) {
      min_row = -1*min_row
      imputed.data[met, is.na(data[met,])] = tryCatch(-1*runif(sum(is.na(data[met,])), min = 0.99*min_row, max= min_row), 
                                                      error = function(e) e, warning=function(w) print(sprintf("%s: met%d", w, met)))
    } else {
      imputed.data[met, is.na(data[met,])] = tryCatch(runif(sum(is.na(data[met,])), min = 0.99*min(rowData), max= min(rowData)), 
                                                      error = function(e) e, warning=function(w) print(sprintf("%s: met%d", w, met)))
    }
  }
  return(imputed.data)
}

zscoreData = function(data, ref) {
  print("zscoreData() called.")
  
  # Only metabolites that also occur in the reference population can be z-scored
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = data[sort(rownames(data)),]
  ref = ref[sort(rownames(ref)),]
  
  # Log transform data
  data = log(data)
  ref = log(data.matrix(ref))
  
  zscore.data = matrix(NA, nrow=nrow(data), ncol=ncol(data))
  rownames(zscore.data) = rownames(data)
  colnames(zscore.data) = colnames(data)
  for (met in 1:nrow(data)) {
    met_data = as.numeric(ref[met,])
    rmSamples = unique(c(which(is.na(met_data)), which(is.infinite(met_data))))
    if (length(rmSamples)>0) {
      x = met_data[-rmSamples]
    } else {
      x = met_data
    }
    if (all(is.na(x))) {
      
    } else {
      #x = x[intersect(which(x>quantile(x, 0.025)), which(x<quantile(x, .975)))]
      d = qqnorm(x, plot.it = FALSE);
      x = as.numeric(d$y)
      z = as.numeric(d$x)
      df = data.frame(x=x,z=z)
      t = lm(x~z, data=df)
      mn.est = as.numeric(t$coefficients[1])
      sd.est = as.numeric(t$coefficients[2])
      rm(d,x,z,df,t)
      zscore.data[met,] = (data[met, ]-mn.est)/sd.est
    }
  }
  
  return(zscore.data)
}

formatNumDataset=function(x){
  x=cleanUpRownames(x,useRownames = F)
  x=apply(x, c(1,2), as.numeric)
}


combineDatasets = function(ref, research) {
  print("combineDatasets() called.")
  ref = as.matrix(ref)
  research = as.matrix(research)
  
  unionMets = unique(c(rownames(ref), rownames(research)))
  data = matrix(0, nrow=length(unionMets), ncol=ncol(ref)+ncol(research))
  rownames(data) = unionMets
  colnames(data) = c(colnames(ref), colnames(research))
  for (r in 1:length(unionMets)) {
    if (unionMets[r] %in% rownames(ref)) {
      data[r,colnames(ref)] = as.numeric(ref[unionMets[r], ])
    } else {
      data[r,colnames(ref)] = rep(NA, ncol(ref))
    }
    
    if (unionMets[r] %in% rownames(research)) {
      data[r,colnames(research)] = as.numeric(research[unionMets[r], ])
    } else {
      data[r,colnames(research)] = rep(NA, ncol(research))
    }
  }
  
  return(data)
}

