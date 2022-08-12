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
#'
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

#' #' plot blown-out modules in graph
#' #'
#' #' @param ig - igraph
#' #' @param S - a set of node names
#' #' @param header4x - column name of sampelAttr used for x-axis categories.
#' #' @param header4color - column name of sampelAttr used for color filling
#' #' @param header4ID -column name of sampelAttr indicating sample IDs.
#' #' @param palette color palette
#' #' @import visNetwork
#' #' @importFrom reshape2 melt
#' #' @example
#' #' @export
#' plotGraph=function(ig,S,size=NULL,membership,color=NULL,visNetworkOP=F,showLable=T){
#'   if()
#'
#'
#' }
