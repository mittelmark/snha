\name{mgraph_components} 
\alias{mgraph_components} 
\title{extract graph components} 
\description{ 
   The function `mgraph_components` extracts the components of a graph 
   which are the clusters of nodes which are not connected to each other. 
} 
\usage{mgraph_components(A)} 
\arguments{ 
\item{A}{either an adjacency matrix of an undirected or direct graph, or an snha object} 
} 
\value{return a vector of component ids for each node} 
\examples{ 
set.seed(123) 
A <- matrix(c(0,1,0, 1,0,0, 0,0,0), nrow=3, byrow=FALSE)  # replace with your adjacency 
mgraph_components(A) 
A = matrix(rbinom(100,1,p=0.1),ncol=10) 
A=A+t(A) ; A[A>1]=1 
mgraph_components(A) 
} 
 
