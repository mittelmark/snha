\name{mgraph} 
\alias{mgraph} 
\title{create a mgraph object, an adjacency matrix of a specific type} 
\description{ 
    This is a function to create a few different graph types such as barabasi, random, cluster etc to 
    test the snha algorithms. 
} 
\usage{mgraph( 
  x=NULL, 
  type="random", 
  mode="directed", 
  nodes=10, 
  edges=12, 
  m=1, 
  k=3, 
  p=NULL, 
  power=1) 
} 
\arguments{ 
\item{x}{either a adjacency matrix or an adjacency list, if given type is not used} 
\item{type}{graph type, one of 'angie', 'band', 'barabasi', 'circle', 'cluster', 
      'hubs', 'random', 'regular' or 'werner', default: 'random'} 
\item{mode}{either 'undirected' or 'directed', default: 'directed'} 
\item{nodes}{number of nodes for a given type, default: 10} 
\item{edges}{number of edges for a given type, not used for type "barabasi", default: 12} 
\item{m}{number of edges added in each iteration of type is 'barabasi', default: 1} 
\item{k}{the degree for a regular graph, or the number of clusters for a cluster graph or the number of hubs for a hubs graph, default: 2} 
\item{p}{the probabilty for edges if type is 'gnp',default=NULL} 
\item{power}{the power for preferential attachment if type is 'barabasi',  
      1 is linear preferential attachment, below one is sublinear, smaller hubs, 0 is no hubs, 
      above 1 is super linear with super hubs, default: 1} 
} 
\keyword{network} 
\keyword{graph} 
\value{mgraph object, which contains an adjacency matrix} 
\examples{ 
M <- matrix(0,nrow=7,ncol=7) 
rownames(M)=colnames(M)=LETTERS[1:7] 
M[c('A','B'),'C']=1 
M['C','D']=1 
M['D',c('E','F')]=1 
M['E','F']=1 
G=mgraph(M) 
set.seed(125) 
R = mgraph(type="random",nodes=8,edges=9) 
summary(R) 
A = mgraph(type="angie",nodes=8,edges=9) 
A 
C = mgraph(type="circle",nodes=8) 
B = mgraph(type="band",nodes=8) 
W=mgraph(type='werner') 
W 
B=mgraph(type='barabasi') 
B 
} 
\seealso{ \link[snha:snha]{snha}, \link[snha:plot.snha]{plot.snha}} 
 
