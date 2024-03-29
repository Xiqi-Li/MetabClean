% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggUMAP.R
\name{plotUmap}
\alias{plotUmap}
\title{plotUmap}
\usage{
plotUmap(
  data_mx,
  sampleAttr = data.frame(),
  label,
  nGroupMax = 10,
  axistextsize = 10,
  dotsize = 2,
  legendtextsize = 9,
  textlabelsize = 0.5
)
}
\arguments{
\item{data_mx}{\itemize{
\item Normalized, imputed, z-scored data. Data matrix includes features as rows, samples as columns.
}}

\item{sampleAttr}{\itemize{
\item a data frame of sample meta data, whose order by row matches colnames of data_mx
}}

\item{label}{\itemize{
\item column name (or index) of "sampleAttr" that represents sample IDs.
}}

\item{nGroupMax}{\itemize{
\item maximum number of categories for categorical meta data.
}}

\item{axistextsize}{\itemize{
\item axis text size
}}

\item{dotsize}{\itemize{
\item dot size
}}

\item{legendtextsize}{\itemize{
\item legend text size
}}

\item{textlabelsize}{\itemize{
\item text label size
}}
}
\value{
learned network in igraph object
}
\description{
ggplot wrapper for UMAP.
If number of features exceed 10000, 10\% of all features were randomly selected then ordered by level of variance. Out of them, top 1000 features were used.
}
