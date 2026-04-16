\name{mgraph_shortest_paths} 
\alias{mgraph_shortest_paths} 
\title{caluclate shortest paths for a graph using BFS search} 
\description{ 
  This function does a shortest path search for the given matrix using sime and slow BFS searrch. 
} 
\usage{mgraph_shortest_paths(g,mode="directed")} 
\arguments{ 
\item{g}{a mgraph object or an adjacency matrix} 
\item{mode}{either 'directed' or 'undirected', default: 'directed'} 
#  } 
\examples{ 
par(mfrow=c(1,2),mai=rep(0.1,4)) 
A=mgraph(type="random",nodes=6,edges=8) 
mgraph_shortest_paths(A) 
} 
 
