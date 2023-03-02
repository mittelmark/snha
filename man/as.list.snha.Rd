\name{as.list.snha} 
\alias{as.list.snha} 
\title{return a list representation for an snha graph object} 
\description{`as.list.snha` provides a S3 method to convert an snha graph object 
             into a list object which can be for instance used to write a report 
              using the library openxlsx.} 
\arguments{ 
\item{x}{ a snha graph object created with the snha function} 
\item{\ldots}{ additional arguments, delegated to the list command} 
} 
\value{list object with the components 'theta', 'sigma', 'p.value'} 
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
\seealso{\link[snha:plot.snha]{plot.snha}} 
 
