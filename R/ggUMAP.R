#### UMAP -GGPLOT2 WRAPPER ####

#' plotUmap
#'
#' ggplot wrapper for UMAP
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param sampleAttr - a data frame of sample meta data, whose order by row matches colnames of data_mx
#' @param nGroupMax - maximum number of categories for categorical meta data.
#' @param label - column name (or index) of "sampleAttr" that represents sample IDs.
#' @param axistextsize - axis text size
#' @param dotsize - dot size
#' @param legendtextsize - legend text size
#' @param textlabelsize - text label size
#' @return learned network in igraph object
#' @import umap ggplot2
#' @export
plotUmap=function(data_mx,sampleAttr=data.frame(),label,nGroupMax=10,axistextsize=10,dotsize=2,legendtextsize=9,textlabelsize=0.5){
  n_neighbors=round(0.25*as.numeric(ncol(data_mx)))
  min_dist=0.1
  config=umap.defaults
  config$n_neighbors=n_neighbors
  config$min_dist=min_dist

  message("***UMAP scoring ***")
  UMAP=umap(t(as.matrix(data_mx)),config = config)
  scores = data.frame(UMAP$layout)

  message("***UMAP wrapper function***")
  umap.plots=list()
  sampleAttr$colSum=colSums(data_mx)
  sampleAttr$max=apply(data_mx,2,max)
  sampleAttr$median=apply(data_mx,2,median)
  feature.names=apply(sampleAttr, c(2), function(x) length(table(x)))
  feature.names=feature.names[class(sampleAttr[,names(feature.names)])== "numeric"|(feature.names<=nGroupMax&feature.names>1)]
  feature.names=names(feature.names)
  for (feature.name in feature.names){

    if(class(sampleAttr[,feature.name])== "numeric") {
      feature = sampleAttr[,feature.name]
    }else{
      feature = as.factor(sampleAttr[,feature.name])
    }
    scores$feature=feature
    umap.plot=
      ggplot(data = scores, aes(x = X1, y = X2, label=sampleAttr[,label])) +
      geom_point(aes(colour = feature), size = dotsize) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize,colour = "black"),
            axis.text.x = element_text(size = axistextsize,colour = "black"),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize),
            legend.title = element_text(size = legendtextsize),
            legend.text = element_text(size = legendtextsize)) +
      labs(colour = feature.name) +
      geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
    if(class(sampleAttr[,feature.name])== "numeric") {
      umap.plot=umap.plot +
        scale_colour_gradient(low = "blue", high = "red")
    }
    umap.plots[[feature.name]]=umap.plot
  }
  return(umap.plots)
}
