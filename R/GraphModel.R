#' learn partial correlation graph
#'
#' learn partial correlation-based graphic models
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param cases - vector of index or column names corresponding to cases
#' @param controls - vector of index or column names corresponding to controls
#' @param fillRateThreashold - percentage of allowed minimum non-NA value percentage per feature.
#' @return learned network in igraph object
#' @importFrom huge huge huge.select
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
learnPartialCorrelationGraph=function(data_mx,cases,controls,fillRateThreashold){
  print("learnPartialCorrelationGraph() called.")
  ref_data=data_mx[,controls]
  ref_data=ref_data[getFill(ref_data)>fillRateThreashold,]
  case_data=data_mx[,cases]
  case_data=case_data[getFill(case_data)>fillRateThreashold,]
  ref_data=ref_data[intersect(rownames(ref_data),rownames(case_data)),]
  case_data=case_data[intersect(rownames(ref_data),rownames(case_data)),]
  data_mx=data.surrogateProfiles(data = case_data, std = 1, ref_data = ref_data)
  data_mx = data_mx[,-which(duplicated(colnames(data_mx)))]
  inv_covmatt = huge(t(data_mx), method="glasso")
  inv_covmatt_select = huge.select(inv_covmatt, criterion = "stars")
  inv_covmat = as.matrix(inv_covmatt_select$icov[[which(inv_covmatt_select$lambda==inv_covmatt_select$opt.lambda)]])
  diag(inv_covmat) = 0;
  colnames(inv_covmat) = rownames(data_mx)
  ig = graph.adjacency(as.matrix(inv_covmat), mode="undirected", weighted=TRUE, add.colnames='name')
  V(ig)$name = rownames(data_mx)
  print(ig)
  return(ig)
}

#' learn condition-specific models
#'
#' learn condition-specific models
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param cases - vector of index or column names corresponding to cases
#' @param controls - vector of index or column names corresponding to controls
#' @param fillRateThreashold - percentage of allowed minimum non-NA value percentage per feature.
#' @param useIgRef - Boolean. If TRUE, provide ig_ref as graph learned from references
#' @param ig_ref - graph learned from references
#' @return learned network in igraph object
#' @importFrom huge huge huge.select
#' @importFrom CTD data.surrogateProfiles graph.naivePruning
#' @examples
#' data(Miller2015)
#' data_mx.og = as.matrix(Miller2015[,grep("IEM_", colnames(Miller2015))])# one sample per column, one metabolite per row.
#' diagnoses=data_mx.og["diagnosis",]
#' cohorts = list()
#' cohorts$msud = colnames(data_mx.og[,which(data_mx.og["diagnosis",]=="Maple syrup urine disease")])
#' cohorts$pku = colnames(data_mx.og[,which(data_mx.og["diagnosis",]=="Phenylketonuria")])
#' cohorts$ref = colnames(data_mx.og[,which(data_mx.og["diagnosis",]=="No biochemical genetic diagnosis")])
#' data_mx=apply(data_mx.og[unlist(cohorts)],c(1,2),as.numeric)
#' ig_ref=learnPartialCorrelationGraph(data_mx, cohorts$ref, cohorts$ref, fillRateThreashold=0.8)
#' ig_PKU=learnModel(data_mx, cohorts$pku, cohorts$ref, fillRateThreashold=0.8, useIgRef=T, ig_ref=ig_ref)
#' @export
learnModel=function(data_mx,cases,controls,fillRateThreashold,useIgRef=F,ig_ref=NULL){
  print("learnModel() called.")
  print("Learning cases vs control network...")
  ig=learnPartialCorrelationGraph(data_mx,controls,controls,fillRateThreashold)
  if(useIgRef==T){
    print("Using provided control-only network...")
    ig_ref=ig_ref
  }else{
    print("Learning control-only network...")
    ig_ref=learnPartialCorrelationGraph(data_mx,cases,controls,fillRateThreashold)
  }
  print("Pruning for condition-specific network...")
  ig_pruned = CTD::graph.naivePruning(ig, ig_ref)
  return(ig_pruned)
}
