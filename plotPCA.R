#' learn partial correlation graph
#'
#' learn partial correlation-based graphic models
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param cases - vector of index or column names corresponding to cases
#' @param controls - vector of index or column names corresponding to controls
#' @param fillRateThreashold - percentage of allowed minimum non-NA value percentage per feature.
#' @return learned network in igraph object
#' @import plotly
#' @importFrom stats prcomp
#' @importFrom CTD data.surrogateProfiles
#' @examples
#' data(Miller2015)
#' data_mx.og = as.matrix(Miller2015[,grep("IEM_", colnames(Miller2015))])# one sample per column, one metabolite per row.
#' diagnoses=data_mx.og["diagnosis",]
#' cohorts = list()
#' cohorts$msud = colnames(data_mx.og[,which(data_mx.og["diagnosis",]=="Maple syrup urine disease")])
#' cohorts$pku = colnames(data_mx.og[,which(data_mx.og["diagnosis",]=="Phenylketonuria")])
#' cohorts$ref = colnames(data_mx.og[,which(data_mx.og["diagnosis",]=="No biochemical genetic diagnosis")])
#' data_mx=apply(data_mx.og[unlist(cohorts)],c(1,2),as.numeric)
#' ig=learnPartialCorrelationGraph(data_mx, cohorts$pku, cohorts$ref, fillRateThreashold=0.8)
#' @export
## PCA
p.pca3=list()
p.pca2=list()
for (tissue in c("AF")){
  load("AFdata.060421.RData")
  for (level in c("zscore","scaledImputed")){
    data_mx=dataSet[[level]]
    data_mx=data_mx[,sampleAttr$PARENT_SAMPLE_NAME]
    features=c("GESTATIONAL_AGE_DAYS","GROUP1","GROUP2")
    temp=prcomp(t(data_mx[getFill(data_mx)==1,]))
    for (feature in features){
      if (feature=="GESTATIONAL_AGE_DAYS"){
        sampleAttr=sampleAttr[sampleAttr$GESTATIONAL_AGE!="UNKNOWN",]
        temp$x=temp$x[sampleAttr$PARENT_SAMPLE_NAME,]
      } else if (feature=="GROUP2"){
        sampleAttr=sampleAttr[sampleAttr$GROUP2 %in% c("normal_M","normal_F","47_XX_21"),]
        temp$x=temp$x[sampleAttr$PARENT_SAMPLE_NAME,]
      }
      p.pca3[[tissue]][[level]][[feature]] = plot_ly(x=temp$x[,1], y=temp$x[,2], z=temp$x[,3], type="scatter3d", mode="markers",
                                                     color = sampleAttr[,feature],text=sampleAttr$PARENT_SAMPLE_NAME) %>%
        layout(title = sprintf("3D PCA of %s Samples by %s (%s)",tissue,capitalize(tolower(feature)),level),scene=list(camera=list(eye=list(x = 1, y =-1, z = 1.5))))
      if (feature=="GESTATIONAL_AGE_DAYS" & level == "zscore") {
        pc=3
      } else {
        pc=2}
      p.pca2[[tissue]][[level]][[feature]] = plot_ly(x=temp$x[,1], y=temp$x[,pc], type="scatter", mode="markers",
                                                     color = sampleAttr[,feature],text=sampleAttr$PARENT_SAMPLE_NAME,
                                                     width = 4) %>%
        layout(title = sprintf("2D PCA of %s Samples by %s",tissue,capitalize(tolower(feature)),level),
               # legend = list(x = 0.7, y = 0.9),
               scene=list(
                 aspectration=list(x=1,y=1)))
      # orca(p,sprintf("%s/PCA.tissue%s.%s.nSample%s.%s.pdf",OP.dir,tissue,feature,dim(data_mx)[2]))
    }
  }
}

#' plot correlation between phenotypes or meta-level traits
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
#'
#' @import reshape2 melt
#'

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
