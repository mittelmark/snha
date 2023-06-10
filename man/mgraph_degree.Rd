\name{mgraph_degree} 
\alias{mgraph_degree} 
\title{return the number of edges for each node} 
\description{ 
    This function returns the number of edges adjacent to every node, 
    for directed graphs it will as well return the in and out going edges. 
} 
\usage{mgraph_degree( 
  g, 
  mode='undirected') 
} 
\arguments{ 
\item{g}{a mgraph object} 
\item{mode}{how should the graph been analyzed, either 'undirected', 'out' or 'in', default: 'undirected'} 
#  } 
\examples{ 
A=mgraph(type="regular",nodes=8,k=3) 
mgraph_degree(A) 
mgraph_degree(A,mode="out") 
mgraph_degree(A,mode="in") 
} 
 
