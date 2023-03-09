\name{as.list.snha} 
\alias{as.list.snha} 
\title{return a list representation for an snha graph object} 
\usage{\method{as.list}{snha}(x,...)} 
\description{The function `as.list.snha` provides a S3 method to convert 
    a snha graph object into a list object which can be for instance used  
    to write a report into an XLSX file  using the library openxlsx.} 
\arguments{ 
\item{x}{ snha graph object created with the snha function} 
\item{\ldots}{ additional arguments, delegated to the list command} 
} 
\value{list object with the components:  
  'chains' (the association chain), 
  'data' (original data), 
  'theta' (adjacency matrix,  
  'sigma' (correlations),  
  'p.value' (correlation p-values) 
} 
\examples{  
data(swiss) 
as=snha(swiss,method="spearman",alpha=0.1) 
result=as.list(as) 
ls(result) 
result$settings 
# can be writte as xlsx file for instance like: 
# library(openxlsx) 
# write.xlsx(result,file="some-result.xlsx") 
} 
\seealso{\link[snha:plot.snha]{plot.snha}, \link[snha:snha]{snha}} 
 
