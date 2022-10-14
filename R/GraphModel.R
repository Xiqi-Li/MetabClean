#' learn partial correlation graph
#'
#' learn partial correlation-based graphic models
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param cases - vector of index or column names corresponding to cases
#' @param controls - vector of index or column names corresponding to controls
#' @param fillRateThreashold - percentage of allowed minimum non-NA value percentage per feature in cases and controls.
#' @return learned network in igraph object
#' @importFrom huge huge huge.select
#' @import igraph
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
  data_mx=CTD::data.surrogateProfiles(data = case_data, std = 1, ref_data = ref_data)
  if(sum(duplicated(colnames(data_mx)))){
    data_mx = data_mx[,-which(duplicated(colnames(data_mx)))]
  }
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
#' @import igraph
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
  ig=learnPartialCorrelationGraph(data_mx,cases,controls,fillRateThreashold)
  if(useIgRef==T){
    print("Using provided control-only network...")
    ig_ref=ig_ref
  }else{
    if(!is.null(ig_ref)){message("User provided reference graph not used; overode by useIgRef")}
    print("Learning control-only network...")
    ig_ref=learnPartialCorrelationGraph(data_mx,controls,controls,fillRateThreashold)
  }
  print("Pruning for condition-specific network...")
  ig_pruned = CTD::graph.naivePruning(ig, ig_ref)
  return(ig_pruned)
}

#' getCTDmodule
#' @example
#' MOI=c(names(orderNCut(data_mx[,cohorts$MSD[1]],threshold=1.7,N=2*kmx)),
#' names(orderNCut(data_mx[,cohorts$MSD[2]],threshold=1.7,N=2*kmx)),
#' names(orderNCut(rowMeans(data_mx[,cohorts$MSD]),threshold=1.7,N=2*kmx))
#' )
#' ranks=getRanks(S=unique(MOI),ig=ig)
#' @export
getRanks=function(S,ig,p1=0.9,thresholdDiff=0.01){
  ## get adjacency_matrix and G
  adjacency_matrix = as.matrix(get.adjacency(ig, attr="weight"))
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name
  S=S[S%in%names(G)]
  ## ranks ##
  ranks = list()
  for (n in 1:length(S)) {
    ind = which(names(G)==S[n])
    #probability_difussion starting from node with index ind
    ranks[[S[n]]]=CTD::singleNode.getNodeRanksN(ind,G,p1=p1,
                                                thresholdDiff=thresholdDiff,
                                                adj_mat=adjacency_matrix,S=S,
                                                num.misses=log2(length(G)),TRUE)
  }
  return(ranks)
}

#' getCTDmodule
#'
#' get best compressed subset of input nodes and out put significance
#' @param ig - background network (condition-specific network).
#' @param profile - named vector of sample omic profile named by corresponding feature name. If useS set as TRUE and S is provided, kmx can be NULL.
#' @param kmx - number of top perturbed features to use as input node set. If useS set as TRUE and S is provided, kmx can be NULL.
#' @param useRanks - Boolean. Whether to use user provided ranks.
#' @param ranks - if useRanks is set as TRUE, provide calculated ranks here.
#' @param useS - Boolean. If TRUE, provide S as input node set.
#' @param S - vector of input node set.
#' @param p1 - parameter passed on to singleNode.getNodeRanksN. if useRanks is TRUE and ranks is provided then this value is ignored.
#' @param thresholdDiff - parameter passed on to singleNode.getNodeRanksN. if useRanks is TRUE and ranks is provided then this value is ignored.
#' @return a list of best compressed node set, p value.
#' @importFrom huge huge huge.select
#' @importFrom CTD data.surrogateProfiles graph.naivePruning mle.getEncodingLength
#' @examples
#' @export
getCTDmodule=function(ig,profile,kmx=30,useRanks=T,ranks=NULL,useS=F,S=NULL,p1=0.9,thresholdDiff=0.01){
  require(igraph)
  ## get S
  if(!useS){
    profile = profile[which(names(profile) %in% V(ig)$name)]
    S = profile[order(abs(profile), decreasing = TRUE)][1:kmx]
    S = names(S)
  }else{
    if(is.null(S)){
      message('Please provide node set S if useS is TRUE.')
    }
  }
  ## get adjacency_matrix and G
  adjacency_matrix = as.matrix(get.adjacency(ig, attr="weight"))
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name

  ## ranks ##
  if(!useRanks){

    message('Pre-computed ranks not provided.
            Getting the single-node encoding node ranks starting from each node...')
    ranks = list()
    for (n in 1:length(S)) {
      ind = which(names(G)==S[n])
      #probability_difussion starting from node with index ind
      ranks[[S[n]]]=CTD::singleNode.getNodeRanksN(ind,G,p1=p1,
                                                  thresholdDiff=thresholdDiff,
                                                  adj_mat=adjacency_matrix,S=S,
                                                  num.misses=log2(length(G)),TRUE)
      # ranks[[S[n]]]=CTD::singleNode.getNodeRanksN(ind,G,p1=p1,
      #                                             thresholdDiff=thresholdDiff,
      #                                             adj_mat=adjacency_matrix,
      #                                             S =NULL,
      #                                             num.misses=NULL,
      #                                             TRUE)

    }
  }else{
    if(is.null(ranks)){
      message('Please provide Pre-computed ranks if useRanks is T.')
    }
  }

  # compute Module and p-value
  S=S[S %in% names(ranks)]
  ptBSbyK = CTD::mle.getPtBSbyK(S, ranks,num.misses = log2(length(G))) # have to specify num.misses to make sure part of ranks work as same as full ranks
  res = CTD::mle.getEncodingLength(ptBSbyK, NULL, NULL, G)
  best_module=ptBSbyK[[which.max(res$d.score)]]
  bestMod=names(best_module)[best_module==1]
  bestCompressedNodeSet=S[S %in% bestMod]
  if(max(res$d.score)<0){max(res$d.score)=0}
  result=list(`compressed node set`=bestCompressedNodeSet,
              `p-value`=2^-(max(res$d.score)),
              `bits`=max(res$d.score))
  if(!useS){
    result[["profile value"]]=profile[bestCompressedNodeSet]
    result[["profile value of input S"]]=profile[S]
    }

  return(result)
}

#' get distance to disease module (or other node set)
#'
#' @param dis_mod - a vector of disease module
#' @example
#' @export
disFromDowntown = function(dis_mod, ptBSbyK.dis, p2.sig.nodes, p2.optBS, ranks, G) {
  p1.e = mle.getEncodingLength(ptBSbyK.dis, NULL, ptID, G)[,"IS.alt"]
  p2.e = mle.getEncodingLength(p2.optBS, NULL, ptID2, G)[,"IS.alt"]

  # Using predefined node ranks, get optimal bitstring for encoding of patient1's union patient2's subsets.
  p12.sig.nodes = unique(c(dis_mod, p2.sig.nodes))
  p12.e = c()
  for (k in 1:length(ptBSbyK.dis)) {
    dis_mod_cpy = dis_mod
    p2.sig.nodes_cpy = p2.sig.nodes

    dis_mod_k = names(which(ptBSbyK.dis[[k]]==1))
    p2.sig.nodes_k = names(which(p2.optBS[[k]]==1))
    while (length(dis_mod_k)<k) {
      dis_mod_k = unique(c(dis_mod_k, dis_mod_cpy[1]))
      dis_mod_cpy = dis_mod_cpy[-1]
    }
    while (length(p2.sig.nodes_k)<k) {
      p2.sig.nodes_k = unique(c(p2.sig.nodes_k, p2.sig.nodes_cpy[1]))
      p2.sig.nodes_cpy = p2.sig.nodes_cpy[-1]
    }
    p12.sig.nodes_k = sapply(unique(c(dis_mod_k, p2.sig.nodes_k)), trimws)
    p12.optBS = mle.getPtBSbyK(p12.sig.nodes_k, ranks, num.misses = log2(length(G)))
    res = mle.getEncodingLength(p12.optBS, NULL, NULL, G)
    p12.e[k] = res[which.max(res[,"d.score"]), "IS.alt"] + log2(choose(length(G), 1))*(length(p12.sig.nodes_k)-which.max(res[,"d.score"]))
  }
  return (list(p1.e=p1.e, p2.e=p2.e, p12.e=p12.e,
               NCD=(p12.e-apply(cbind(p1.e, p2.e), 1, min))/apply(cbind(p1.e, p2.e), 1, max)))
}

#' get disease module
#'
#' calculate disease module
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param cases - A Vector of index or column names corresponding to cases.
#' @param kmx - The number of top perturbed features to use as input node set. If useS set as TRUE and S is provided, kmx can be NULL.
#' @param zThreshold - z-score threshold where values outside of range (-zThreshold,zThreshold) are considered abnormal.
#' @param thresholdDiff - parameter passed on to singleNode.getNodeRanksN.
#' @param ranksList - A list of ranks. If CrossValidated is TRUE, provide ranks calculated from network folds with names corresponding to the left-out sample, and a ranks object calculatd from the network learned from all training samples named by "0".
#' @param igList - A list of background networks. If CrossValidated is TRUE, provide network folds with names corresponding to the left-out sample, and a network learned from all training samples named by "0".
#' @param CrossValidated - set as TRUE if cross-validated method were used when learning background networks.
#' @param useCasesMean - Use mean cases profile for calculating initial disease module.
#' @return a list of best compressed node set, p value.
#' @importFrom huge huge huge.select
#' @importFrom CTD data.surrogateProfiles graph.naivePruning mle.getPtBSbyK mle.getEncodingLength
#' @examples
#' @export

getDiseaseModule=function(data_mx,cases,kmx=30,zThreshold,ranksList,igList,CrossValidated=F,useCasesMean=T){
  # set function

  data_mx=data_mx[,which(colnames(data_mx) %in% cases)]
  downtown_disease_module = c()

  if(useCasesMean){
    message("Calculating best compressed set from average disease profile")
    mn = apply(data_mx, 1, function(i) mean(na.omit(i)))
    if(sum(is.na(mn))>0){mn = mn[-which(is.na(mn))]}

    if(CrossValidated){
      for (pt in cases) {
        ig = igList[[pt]]
        ranks = ranksList[[pt]]
        adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))

        G = vector(mode="list", length=length(V(ig)$name))
        names(G) = V(ig)$name
        S = mn[which(abs(mn)>zThreshold)] # change mean z-score threshold for abnomral molecules
        S = S[which(names(S) %in% names(G))]
        it=0
        while (length(S)< kmx) { # change the length of wanted abnomral molecules
          S = mn[which(abs(mn)>(zThreshold-0.1*it))]
          S = S[which(names(S) %in% names(G))]
          it = it + 1
        }
        ptBSbyK = mle.getPtBSbyK(names(S), ranks, num.misses = log2(length(G)))
        res = mle.getEncodingLength(ptBSbyK, NULL, NULL, G)
        downtown_disease_module = c(downtown_disease_module, names(which(ptBSbyK[[which.max(res[,"d.score"])]]==1)))
      }
    }else{
      ig = igList[["0"]]
      ranks = ranksList[["0"]]
      adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))

      G = vector(mode="list", length=length(V(ig)$name))
      names(G) = V(ig)$name
      S = mn[which(abs(mn)>zThreshold)] # change mean z-score threshold for abnomral molecules
      S = S[which(names(S) %in% names(G))]
      it=0
      while (length(S)< kmx) { # change the length of wanted abnomral molecules
        S = mn[which(abs(mn)>(zThreshold-0.1*it))]
        S = S[which(names(S) %in% names(G))]
        it = it + 1
      }
      ptBSbyK = mle.getPtBSbyK(names(S), ranks, num.misses = log2(length(G)))
      res = mle.getEncodingLength(ptBSbyK, NULL, NULL, G)
      downtown_disease_module = c(downtown_disease_module, names(which(ptBSbyK[[which.max(res[,"d.score"])]]==1)))
    }
    print(sprintf("Initial size of downtown disease module (calculated from averaged cases profile): %s",length(unique(downtown_disease_module))))
  }else{
    message("Calculating best compressed set from all cases profiles")
    for (pt in cases) {
      if(CrossValidated){
        ig=igList[[pt]]
        ranks=ranksList[[pt]]
      }else{
        ig=igList[["0"]]
        ranks=ranksList[["0"]]
      }

      S=orderNCut(data_mx[rownames(data_mx) %in% V(ig)$name,pt],zThreshold,kmx)
      tmp=data_mx[names(S),pt]
      it=0
      while (length(S)< kmx) { # change the length of wanted abnomral molecules
        S = tmp[which(abs(tmp)>(zThreshold-0.1*it))]
        it = it + 1
      }
      CTD.ind=getCTDmodule(ig=ig,profile=NULL,kmx=NULL,useRanks=T,ranks=ranks,useS=T,S=names(S),p1=0.9,thresholdDiff=0.01)
      downtown_disease_module = c(downtown_disease_module, CTD.ind$`compressed node set`)
    }
  }

  # Estimate disease module "purity" based on 2 size thresholds and known case mean distances
  ig = igList[["0"]]
  ranks = ranksList[["0"]]
  adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))
  G = vector(mode="list", length=length(V(ig)$name))
  names(G) = V(ig)$name
  data_mx = data_mx[which(rownames(data_mx) %in% V(ig)$name), ]
  purity_mets = list()
  purity = c()
  # Size1: Full downtown disease module
  tmp.mm = unique(downtown_disease_module)
  tmp.mm = tmp.mm[which(tmp.mm %in% names(G))]
  ptBSbyK.dis = mle.getPtBSbyK(tmp.mm, ranks, num.misses = log2(length(G)))
  res = mle.getEncodingLength(ptBSbyK.dis, NULL, ptID, G)
  downtown_disease_mod = names(which(ptBSbyK.dis[[which.max(res[,"d.score"])]]==1))
  dd = c()
  for (pt in 1:length(cases)) {
    ptBSbyK.dis = mle.getPtBSbyK(unique(downtown_disease_mod), ranks, num.misses = log2(length(G)))
    p2.sig.nodes = names(data_mx[order(abs(data_mx[,cases[pt]]), decreasing = TRUE), cases[pt]][1:length(unique(downtown_disease_mod))])
    p2.optBS = mle.getPtBSbyK(p2.sig.nodes, ranks, num.misses = log2(length(G)))
    ctdsim = disFromDowntown(unique(downtown_disease_mod), ptBSbyK.dis, p2.sig.nodes, p2.optBS, ranks, G)
    dd[pt] = min(ctdsim$NCD)
  }
  purity_mets[[1]] = unique(downtown_disease_mod)
  purity[1] = mean(dd)
  print(sprintf("Full size of the best compressed downtown disease module: %s; purity: %s",length(unique(downtown_disease_mod)),purity[1]))

  if (!(CrossValidated==F & useCasesMean==T)){
    # Size2: Smaller disease module
    downtown_disease_module = names(which(table(downtown_disease_module)>(length(cohorts[[model]])/2)))
    tmp.mm = unique(downtown_disease_module)
    tmp.mm = tmp.mm[which(tmp.mm %in% names(G))]
    ptBSbyK.dis = mle.getPtBSbyK(unique(downtown_disease_module), ranks, num.misses = log2(length(G)))
    res = mle.getEncodingLength(ptBSbyK.dis, NULL, ptID, G)
    downtown_disease_mod = names(which(ptBSbyK.dis[[which.max(res[,"d.score"])]]==1))
    print(length(unique(downtown_disease_mod)))
    dd = c()
    for (pt in 1:length(cases)) {
      ptBSbyK.dis = mle.getPtBSbyK(unique(downtown_disease_mod), ranks, num.misses = log2(length(G)))
      p2.sig.nodes = names(data_mx[order(abs(data_mx[,pt]), decreasing = TRUE), pt][1:length(unique(downtown_disease_mod))])
      p2.optBS = mle.getPtBSbyK(p2.sig.nodes, ranks, num.misses = log2(length(G)))
      ctdsim = disFromDowntown(unique(downtown_disease_mod), ptBSbyK.dis, p2.sig.nodes, p2.optBS, ranks, G)
      dd[pt] = min(ctdsim$NCD)
    }
    purity_mets[[2]] = downtown_disease_mod
    purity[2] = mean(dd)
    print(sprintf("Smaller size of the best compressed downtown disease module: %s; purity: %s",length(unique(downtown_disease_mod)),purity[2]))
  }


  # Selected disease module based on "purity" estimates
  downtown_disease_module = purity_mets[[which.min(purity)]]
  result=list(disease_module=downtown_disease_module,
              cases_module_profiles=data_mx[downtown_disease_module,],
              features=purity_mets,
              purity=purity)
  return(result)

}

#' get CTD distance to disease module
#' @param data_mx - Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
#' @param cases - A Vector of column names corresponding to cases.
#' @param diseaseModule - A Vector of node names representing disease module.
#' @param igList - A list of graphs with names corresponding to names of Folds as character.
#' @param rankList - A list of ranks with names corresponding to names of Folds as character.
#' @param Folds - A list of assignment of left out samples of each training fold.
#' @import CTD igraph
#' @export
CTDdm=function(data_mx,cases,diseaseModule,igList,rankList,Folds){
  mn = apply(data_mx[,which(colnames(data_mx) %in% cases)], 1, function(i) mean(na.omit(i)))
  if(sum(is.na(mn))>0){mn = mn[-which(is.na(mn))]}

  df_DISMOD=list()
  for (fold in c(0,seq(Folds))){
    ig = igList[[as.character(fold)]]
    ranks = ranksList[[as.character(fold)]]
    adjacency_matrix = list(as.matrix(get.adjacency(ig, attr="weight")))
    G = vector(mode="list", length=length(V(ig)$name))
    names(G) = V(ig)$name
    data_mx = data_mx[which(rownames(data_mx) %in% V(ig)$name), ]

    tmp.mm = unique(diseaseModule)
    tmp.mm = tmp.mm[which(tmp.mm %in% names(G))]
    ptBSbyK.dis = mle.getPtBSbyK(tmp.mm, ranks, num.misses = log2(length(G)))
    res = mle.getEncodingLength(ptBSbyK.dis, NULL, NULL, G)
    downtown_disease_mod = names(which(ptBSbyK.dis[[which.max(res[,"d.score"])]]==1))
    ptBSbyK.dis = mle.getPtBSbyK(downtown_disease_mod, ranks,num.misses = log2(length(G)))

    df_disMod = data.frame(stringsAsFactors = FALSE)
    for (p in 1:ncol(data_mx)) {
      print(sprintf("Patient %d/%d...", p, ncol(data_mx)))
      ptID = colnames(data_mx)[p]
      if (ptID %in% cases) {diag = model}else{diag = "diagnose_not_match_model"}
      # CTD: get module best explained by network's connectivity
      S = data_mx[order(abs(data_mx[,p]), decreasing = TRUE),p]
      p2.sig.nodes = names(S)[1:length(downtown_disease_mod)]
      p2.optBS = mle.getPtBSbyK(p2.sig.nodes, ranks, num.misses = log2(length(G)))
      ctdDisMod = disFromDowntown(downtown_disease_mod, ptBSbyK.dis, p2.sig.nodes, p2.optBS, ranks, G)
      #if (ptID=="X606789") {print(min(ctdDisMod$NCD))}
      df_disMod[p, "ptID"] = colnames(data_mx)[p]
      df_disMod[p, "diag"] = diag[1]
      df_disMod[p, "ctdDisMod"] = min(ctdDisMod$NCD)
    }
    df_DISMOD[[as.character(fold)]]=df_disMod
  }
  return(df_DISMOD)
}




