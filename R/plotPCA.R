#' PCA plot wrapper
#'
#' plot PCA in plotly and ggplot
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param sampleAttr - a data frame of sample meta data, whose order by row matches colnames of data_mx
#' @param nGroupMax - maximum number of categories for categorical meta data.
#' @param label - column name (or index) of "sampleAttr" that represents sample IDs.
#' @param axistextsize - axis text size
#' @param dotsize - dot size
#' @param legendtextsize - legend text size
#' @param textlabelsize - text label size
#' @return list of plotly and/or ggplot object(s)
#' @import plotly ggplot dplyr
#' @importFrom stats prcomp
#'
#' @export
## PCA
plotPCA=function(data_mx,sampleAttr=data.frame(),
                 label,nGroupMax=10,dimensions=c(2,3),
                 axistextsize=10,dotsize=2,legendtextsize=9,textlabelsize=0.5){
  if(nrow(data_mx)>10000){
    message("*** selecting top 1000 variable features ***")
    data_mx=data_mx[
      names(head(sort(apply(data_mx[sample(1:nrow(data_mx),nrow(data_mx)/10),], 1, var),decreasing=T),1000)),]
  }

  #calculate PCA
  temp=prcomp(t(data_mx[getFill(data_mx)==1,]))

  # determine features
  sampleAttr$colSum=colSums(data_mx)
  sampleAttr$max=apply(data_mx,2,max)
  sampleAttr$median=apply(data_mx,2,median)
  feature.names=apply(sampleAttr, c(2), function(x) length(table(x)))
  feature.names=feature.names[sapply(sampleAttr[names(feature.names)],class)== "numeric"|(feature.names<=nGroupMax&feature.names>1)]
  feature.names=names(feature.names)

  # plot pca with color indicative of feature meta info
  p.pca=list()
  for (feature.name in feature.names){
    if(class(sampleAttr[,feature.name])== "numeric") {
      feature = sampleAttr[,feature.name]
    }else{
      feature = as.factor(sampleAttr[,feature.name])
    }
    if(2 %in% dimensions){
      p.pca[["3D"]][[feature.name]] = plot_ly(x=temp$x[,1], y=temp$x[,2], z=temp$x[,3],
                                              type="scatter3d", mode="markers",
                                              color = feature,
                                              # marker = list(
                                              # colorscale='Hot',
                                              # colorbar=list(title='Colorbar')
                                              # ),
                                              text=sampleAttr[,label]) %>%
        layout(title = sprintf("3D PCA by %s",capitalize(tolower(feature.name))),
               scene=list(camera=list(eye=list(x = 1, y =-1, z = 1.5))))
    }
    if(3 %in% dimensions){
      p.pca[["2D"]][[feature.name]] = ggplot(data = temp$x, aes(x = PC1, y = PC2, label=sampleAttr[,label])) +
        geom_point(aes(colour = feature)) +
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
        ggtitle(sprintf("2D PCA by %s",capitalize(tolower(feature.name))))+
        geom_text(vjust = "inward", hjust = "inward", size = textlabelsize)
    }
  return(p.pca)
  }
}


