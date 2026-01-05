\name{mgraph_harmonic_centrality} 
\alias{mgraph_harmonic_centrality} 
\title{return the harmonic centraility for each node of a graph} 
\description{ 
    This function returns the harmonic centrality for every node, of a graph 
    directed graphs are converted to unidirecred graphs before and then the following 
    formula is applied: \code{h = sum(1/d) / (n-1)} where d are the shortest paths 
    of a node to all other nodes of the graph and n ist the number of nodes in a graph. 
} 
\usage{mgraph_harmonic_centrality( 
  g, 
  norm=TRUE) 
} 
\arguments{ 
\item{g}{a mgraph object} 
\item{norm}{should the harmonic centrality normalized over the number of nodes, default: TRUE} 
#  } 
\examples{ 
 A=mgraph(type="regular",nodes=8,k=3) 
 plot(A,layout="sam") 
 mgraph_degree(A) 
 mgraph_harmonic_centrality(A) 
 B=mgraph(type="angie",nodes=20,edges=30) 
 plot(B,layout="sam") 
 mgraph_degree(B) 
 mgraph_harmonic_centrality(B) 
} 
 
