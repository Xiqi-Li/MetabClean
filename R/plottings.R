#' plot correlation between phenotypes or meta-level traits
#'
#' @importFrom Hmisc rcorr
#' @param sampleAttr metadata table with phenotype data listed by columns
#' @param type correlation method
#' @importFrom ggcorrplot ggcorrplot
#' @export
plotCor=function(sampleAttr,type=c("pearson","spearman")){
  cormx=Hmisc::rcorr(apply(sampleAttr,c(1,2),as.numeric),type=c("pearson"))
  p=ggcorrplot::ggcorrplot(cormx$r,hc.order = TRUE,type = "lower",lab = TRUE,method = c("circle"),p.mat = cormx$P)
  return(p)
}

#' plot box plots
#'
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param sampleAttr - sample attributes; meta data
#' @param MOI - molecule of interest. A subset of data_mx colnames.
#' @param header4x - column name of sampelAttr used for x-axis categories.
#' @param header4color - column name of sampelAttr used for color filling
#' @param header4ID -column name of sampelAttr indicating sample IDs.
#' @param palette color palette
#' @importFrom ggpubr ggboxplot
#' @export
plotBoxPlot=function(data_mx,sampleAttr,MOI,header4x,header4color=NULL,header4ID,palette= "Dark2"){
  df=as.data.frame(t(data_mx[MOI,]))
  df[,header4ID]=rownames(df)
  df[,header4x]=sampleAttr[match(df[,header4ID],sampleAttr[,header4ID]),header4x]
  if(!is.null(header4color)){
    df[,header4color]=sampleAttr[match(df[,header4ID],sampleAttr[,header4ID]),header4color]
    df=df[!is.na(df[,header4x])&!is.na(df[,header4color]),]
    p=ggpubr::ggboxplot(df, x = header4x ,
                        y = MOI,
                        combine = TRUE,
                        # merge = "flip",
                        ylab = "Value",
                        add = "jitter",                               # Add jittered points
                        add.params = list(size = 0.5, jitter = 0.1),
                        color = header4color ,
                        palette = "Dark2")
  }else{
    df=df[!is.na(df[,header4x]),]
    p=ggpubr::ggboxplot(df, x = header4x ,
                        y = MOI,
                        combine = TRUE,
                        # merge = "flip",
                        ylab = "Value",
                        add = "jitter",                               # Add jittered points
                        add.params = list(size = 0.5, jitter = 0.1),
                        palette = "Dark2")
  }
  return(p)
}

#' plot scatter plots
#'
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param sampleAttr - sample attributes; meta data
#' @param MOI - molecule of interest. A subset of data_mx colnames.
#' @param header4x - column name of sampelAttr used for x-axis categories.
#' @param header4color - column name of sampelAttr used for color filling
#' @param header4ID -column name of sampelAttr indicating sample IDs.
#' @param palette color palette
#' @importFrom reshape2 melt
#' @example
#' @export

plotScatterPlot=function(data_mx,sampleAttr,MOI,header4x,header4color=NULL,header4ID,palette="Dark2"){
  df=as.data.frame(t(data_mx[MOI,]))
  df[,header4ID]=rownames(df)
  df[,header4x]=sampleAttr[match(df[,header4ID],sampleAttr[,header4ID]),header4x]
  if(!is.null(header4color)){
    df[,header4color]=sampleAttr[match(df[,header4ID],sampleAttr[,header4ID]),header4color]
    df = reshape2::melt(df,id.vars=c(header4x,header4ID,header4color))
    if(class(sampleAttr[,header4color])%in%c("character","factor")){
      p = ggpubr::ggscatter(df, x = header4x, y = "value", color = header4color,
                            palette = palette,
                            add = "reg.line",
                            add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                            conf.int = TRUE, # Add confidence interval
                            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                            facet.by="variable")
    }else{
      p = ggpubr::ggscatter(df, x = header4x, y = "value", color = header4color,
                            # palette = palette,
                            add = "reg.line",  # Add regressin line
                            add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                            conf.int = TRUE, # Add confidence interval
                            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                            cor.coeff.args = list(method = "pearson", label.x = 0, label.sep = "\n"),
                            facet.by="variable")
    }
  }else{
    df = reshape2::melt(df,id.vars=c(header4x,header4ID))
    p = ggpubr::ggscatter(df, x = header4x, y = "value",
                          # palette = palette,
                          add = "reg.line",  # Add regressin line
                          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE, # Add confidence interval
                          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                          cor.coeff.args = list(method = "pearson", label.x = 0, label.sep = "\n"),
                          facet.by="variable")

  }
  p=p+
    theme_bw()+
    theme(axis.text.x=element_blank(),legend.position = "bottom")+
    xlab(header4x)+ylab("value")
  return(p)
}

# plot blown-out modules in graph
#
# @param ig - igraph
# @param S - a set of node names
# @param header4x - column name of sampelAttr used for x-axis categories.
# @param header4color - column name of sampelAttr used for color filling
# @param header4ID -column name of sampelAttr indicating sample IDs.
# @param palette color palette
# @import visNetwork
# @importFrom reshape2 melt
# @example
plotGraph=function(ig,S,size=NULL,membership,color=NULL,visNetworkOP=F,showLable=T){
}





#' plot heat map
#'
#' @param data_mx - data matrix
#' @param classFactor - a factor of categories designation of samples (or column names of data_mx)
#' @param classColorMatch - a named vector of colors, each named by one of all possible categories in classFactor
#' @param colorRowSide - Boolean. Whether heat map row side should be colored.
#' @param rowClassFactor - a factor of categories designation of features (or row names of data_mx)
#' @param rowClassColorMatch - a named vector of colors, each named by one of all possible categories in rowClassFactor
#' @param heatColors - color of heat.
#' @param colorby - heat should by colored by values or ranks of values.
#' @param colorThresh - if colorby == "values", provide values of color change thresholds.
#' @param dendrogram - same as heatmap.2. Show dendrogram on "both", "row", "column" or "none".
#' @param Rowv - same as heatmap.2. cluster row.
#' @param Colv -same as heatmap.2. cluster col.
#' @param labRow -same as heatmap.2. show row labels.
#' @param labCol -same as heatmap.2. show column labels.
#' @example
#' @import gplots
#' @import ggplot2
#' @importFrom cols4all c4a
#' @importFrom cowplot get_legend
#' @importFrom grDevices recordPlot
#' @export
plotHeatmap=function(data_mx,
                     classFactor,
                     classColorMatch=NULL,
                     colorRowSide=F,
                     rowClassFactor=NULL,
                     rowClassColorMatch=NULL,
                     heatColors=c("white","steelblue"),
                     colorby=c("value","rank"),
                     colorThresh=NULL,
                     dendrogram = c("both"),
                     Rowv=T,Colv=T,
                     labRow=T,labCol=T
                     ){
  # set column side label colors
  if(is.null(classColorMatch)){
    classColorMatch=cols4all::c4a("wright25",nlevels(classFactor))
    names(classColorMatch)=levels(classFactor)
  }
  ColSideColors=classColorMatch[classFactor]

  # set row side label colors
  if(colorRowSide){
    if(is.null(rowClassFactor)){
      RowSideColors=ColSideColors
      rowClassColorMatch=classColorMatch
    }else{
      if(is.null(rowClassColorMatch)){
        rowClassColorMatch=cols4all::c4a("dark24",nlevels(rowClassFactor))
        names(rowClassColorMatch)=levels(rowClassFactor)
        RowSideColors=rowClassColorMatch[rowClassFactor]
        }
      }
    }else{
      RowSideColors=NULL
    }

  # set heat map colors
  if(length(heatColors)==3){
    message("Using user provided heatCol_breaks")

    if(colorby=="value"){
      if(is.null(colorThresh)){
        message("please provid colorThresh. try colorThresh=c(-2,2)")
      }else{
        col_breaks = c(seq(min(data_mx),colorThresh[1]-0.0001,length=10), # for low
                       seq(colorThresh[1],colorThresh[2]-0.0001,length=10), # for medium
                       seq(colorThresh[2],max(data_mx),length=10)) # for high
      }
    }else if (colorby=="rank"){
      temp=sort(data_mx)
      tmp=floor(length(temp)/3)
      temp=temp[c(1,tmp,2*tmp,length(temp))]
      col_breaks = c(seq(temp[1],temp[2]-0.0001,length=10), # for low
                     seq(temp[2],temp[3]-0.0001,length=10), # for medium
                     seq(temp[3],temp[4],length=20)) # for high
    }
  }else{
    col_breaks=seq(0,1,0.1)
  }
  heatColorGradient=colorRampPalette(heatColors)(n=length(col_breaks)-1)

  # get legend
  colordf = data.frame(x = 1:length(classColorMatch),
                       y = 1:length(classColorMatch),
                       group = names(classColorMatch))
  ggp = ggplot(colordf, aes(x, y, color = group)) +
    geom_point(shape=15) +
    scale_color_manual(values=classColorMatch) +
    theme_bw() +
    theme(legend.title = element_blank())
  colorLegend <- cowplot::get_legend(ggp)

  colordf = data.frame(x = 1:length(rowClassColorMatch),
                       y = 1:length(rowClassColorMatch),
                       group = names(rowClassColorMatch))
  ggp = ggplot(colordf, aes(x, y, color = group)) +
    geom_point(shape=15) +
    scale_color_manual(values=rowClassColorMatch) +
    theme_bw() +
    theme(legend.title = element_blank())
  rowColorLegend <- cowplot::get_legend(ggp)

  # plot
  gplots::heatmap.2(data_mx,
                    dendrogram=dendrogram,
                    Rowv=Rowv,
                    Colv=Colv,
                    trace="none",
                    col=heatColorGradient,
                    breaks=col_breaks,
                    margins=c(10,10),
                    ColSideColors = ColSideColors,
                    RowSideColors = RowSideColors,
                    labRow = labRow,
                    labCol = labCol
                    )
  # gridGraphics::grid.echo()
  # heatmap = grid::grid.grab()
  # p.grid=gridExtra::grid.arrange(heatmap,colorLegend,nrow=1,widths=c(3.5,0.5))
  heatmap=grDevices::recordPlot()
  p.list=list(heatmap=heatmap,legend=list(colorLegend,rowColorLegend))
  return(p.list)
}
