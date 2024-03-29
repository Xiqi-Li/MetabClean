% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GraphModel.R
\name{getDiseaseModule}
\alias{getDiseaseModule}
\title{get disease module}
\usage{
getDiseaseModule(
  data_mx,
  cases,
  kmx = 30,
  zThreshold,
  ranksList,
  igList,
  CrossValidated = F,
  useCasesMean = T
)
}
\arguments{
\item{data_mx}{\itemize{
\item Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
}}

\item{cases}{\itemize{
\item A Vector of index or column names corresponding to cases.
}}

\item{kmx}{\itemize{
\item The number of top perturbed features to use as input node set. If useS set as TRUE and S is provided, kmx can be NULL.
}}

\item{zThreshold}{\itemize{
\item z-score threshold where values outside of range (-zThreshold,zThreshold) are considered abnormal.
}}

\item{ranksList}{\itemize{
\item A list of ranks. If CrossValidated is TRUE, provide ranks calculated from network folds with names corresponding to the left-out sample, and a ranks object calculatd from the network learned from all training samples named by "0".
}}

\item{igList}{\itemize{
\item A list of background networks. If CrossValidated is TRUE, provide network folds with names corresponding to the left-out sample, and a network learned from all training samples named by "0".
}}

\item{CrossValidated}{\itemize{
\item set as TRUE if cross-validated method were used when learning background networks.
}}

\item{useCasesMean}{\itemize{
\item Use mean cases profile for calculating initial disease module.
}}

\item{thresholdDiff}{\itemize{
\item parameter passed on to singleNode.getNodeRanksN.
}}
}
\value{
a list of best compressed node set, p value.
}
\description{
calculate disease module
}
