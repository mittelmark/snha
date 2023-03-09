\name{snha_layout} 
 
\alias{snha_layout} 
\title{Determine graph layouts} 
\usage{snha_layout( 
   A, 
   mode='sam', 
   method='pearson',  
   noise=FALSE,  
   star.center=NULL, 
   interactive=FALSE) 
} 
\description{This function returns xy coordinates for a given input 
      adjacency matrix or snha graph. It is useful if you like to plot 
      the same set of nodes with different edge connections  
      in the same layout.} 
\arguments{ 
\item{A}{an adjacency matrix or an snha graph object} 
\item{mode}{character string for the layout type, can be either  
     'mds' (mds on graph using shortest paths), 'mdsd' (mds on data) 
     'sam' (sammon on graph), 'samd' (sammon on data), 
     'circle', 'grid' or 'star', default: 'sam'} 
\item{method}{method for calculating correlation distance if mode is either 
      'mdsd' or 'samd', default: 'pearson'} 
\item{noise}{should some noise be added, default: FALSE} 
\item{star.center}{the centered node if layout is 'star', must be 
          a character string for the node name, default: NULL} 
\item{interactive}{switch into interactive mode where you can click in the 
       graph and move nodes with two clicks, first selecting the node, second click 
       gives the new coordinates for the node, default: FALSE} 
} 
\keyword{network} 
\keyword{layout} 
\value{matrix with x and y columns for the layout} 
\examples{ 
data(swiss) 
sw.s=snha(swiss,method='spearman') 
sw.p=snha(swiss,method='pearson') 
lay=snha_layout(sw.s,mode='sam') 
plot(sw.s,layout=lay) 
plot(sw.p,layout=lay) 
plot(sw.s,layout='star',star.center='Education') 
rn1=rnorm(nrow(swiss)) 
nswiss=cbind(swiss,Rn1=rn1) 
plot(snha(nswiss,method='spearman'),layout='sam') 
plot(snha(nswiss,method='spearman'),layout='samd', 
  vertex.size=2,vertex.color='beige') 
} 
