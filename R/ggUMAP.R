#### UMAP -GGPLOT2 WRAPPER ####
require(umap)
require(ggplot2)


plotUmap=function(data_mx,sampleAttr=data.frame(),nGroupMax=10,axistextsize=10,dotsize=2,legendtextsize=9,textlablesize=0.5){
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
  feature.names=feature.names[class(sampleAttr[,feature.names])== "numeric"|(feature.names<=nGroupMax&feature.names>1)]
  feature.names=names(feature.names)
  for (feature.name in feature.names){
    
    if(class(sampleAttr[,feature.name])== "numeric") {
      feature = sampleAttr[,feature.name]
    }else{
      feature = as.factor(sampleAttr[,feature.name])
    }
    scores$feature=feature
    umap.plot=
      ggplot(data = scores, aes(x = X1, y = X2, label=sampleAttr$individualID)) +
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
      geom_text(vjust = "inward", hjust = "inward", size = textlablesize)
    if(class(sampleAttr[,feature.name])== "numeric") {
      umap.plot=umap.plot +
        scale_colour_gradient(low = "blue", high = "red")
    }
    umap.plots[[feature.name]]=umap.plot
  }
  return(umap.plots)
}
