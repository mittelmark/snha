\name{mgraph_d2u} 
\alias{mgraph_d2u} 
\title{create an undirected graph out of a directed graph} 
\description{ 
  This function gets an directed graph and convertes all edges from directed ones to undirected ones. 
  The number of edges should stay the same, the edge sign (+ or -) stays the same. 
} 
\usage{mgraph_d2u(g)} 
\arguments{ 
\item{g}{a mgraph object or an adjacency matrix} 
#  } 
\examples{ 
A=mgraph(type="angie",nodes=4,edges=4) 
A 
U=mgraph_d2u(A) 
U 
} 
 
