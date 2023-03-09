\name{snha_rsquare} 
\alias{snha_rsquare} 
\title{linear model based r-square values for given data and graph} 
\usage{snha_rsquare(data,graph=NULL)} 
\description{ The function `snha_rsquare` calculates for given data and  
a graph the covered r-squared values by  a linear model for each node.  
The linear model predicts each node by an additive mode  
using it's neighbor nodes in the graph.} 
\arguments{ 
\item{data}{data matrix or data frame where variables 
are in columns and samples in rows or a snha graph} 
\item{graph}{ graph object or adjacency matrix of an (un)directed graph, not 
needed if data is a snha graph, default: NULL.} 
} 
\value{vector of rsquare values for each node of the graph} 
\examples{  
# random adjacency matrix 
A=matrix(rbinom(100,1, 0.2),nrow=10,ncol=10) 
diag(A)=0 
colnames(A)=rownames(A)=LETTERS[1:10] 
# random data 
data=matrix(rnorm(1000),ncol=10) 
colnames(data)=colnames(A) 
snha_rsquare(data,A) 
# real data 
data(swiss) 
sw.s=snha(swiss,method='spearman') 
rsqs=snha_rsquare(sw.s) 
plot(sw.s,main=paste("r =",round(mean(rsqs,2))), 
   layout='star',star.center='Examination') 
# some colors for r-square values 
vcols=paste("grey",seq(80,40,by=-10),sep="") 
scols=as.character(cut(snha_rsquare(swiss,sw.s$theta), 
   breaks=c(0,0.1,0.3,0.5,0.7,1),labels=vcols)) 
plot(sw.s,main=paste("r =",round(mean(snha_rsquare(swiss,sw.s$theta)),2)), 
   vertex.color=scols ,layout='star',star.center='Examination', 
   vertex.size=10,edge.color=c('black','red'),edge.width=3) 
} 
 
