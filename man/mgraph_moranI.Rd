\name{mgraph_moranI} 
\alias{mgraph_moranI} 
\alias{mgraph_moranI} 
\title{Measure of global spatial autocorrelation adapted for networks} 
\description{ 
  The function mgraph_moranI calculates the global Moran's I for the clustering of node 
  properties, positive values indicate a tendency of neighborhood 
  of nodes where they have similar property values, whereas negative values indicate a tendency of  
  agglomeration of node with high values towards nodes with low vales. 
  values around zero indicate no association between the node values  
  of neighbor nodes.  
} 
\usage{mgraph_moranI(A, x,B=1000)} 
\arguments{ 
\item{A}{adjacency matrix of an undirected graph} 
\item{x}{node values for the nodes of the adjacency matrix} 
\item{B}{number of boostrap samplings to determin the p-values for the Moran I value, default: 1000} 
} 
\value{List with following components: 
  \item{estimate}{gives Moran's I value which should be between -1 and 1} 
  \item{expected}{gives the expected Moran's I value under the assumption that H0 is true} 
  \item{p.value.greater}{proportion of boostrapped moran I values which where larger than the estimate value} 
  \item{p.value.less}{proportion of boostrapped moran I values which where smaller than the estimate value} 
  \item{p.value.two.sided}{two sided p-value calculated on the boostrapped p-values} 
} 
\references{ 
  \itemize{ 
    \item Moran, P. A. (1950). Notes on continuous stochastic phenomena. \emph{Biometrika, 37(1/2), 17-23}. 
  } 
} 
