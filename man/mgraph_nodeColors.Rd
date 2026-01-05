\name{mgraph_nodeColors} 
\alias{mgraph_nodeColors} 
\title{create node colors for directed graphs} 
\description{ 
  This function simplifies automatic color coding of nodes for directed graphs. 
  Nodes will be colored based on their degree properties, based 
  on their incoming and outcoming edges. 
} 
\usage{mgraph_nodeColors(g, 
    col=c("skyblue","grey80","salmon"), 
    type="inout", 
    legend=FALSE 
)} 
\arguments{ 
\item{g}{a mgraph object or an adjacency matrix or an adjacency list} 
\item{col}{default colors for nodes with only incoming, in- and outgoing and only outgoing edges, default: c("skyblue","grey80","salmon")} 
\item{type}{type of color selection, either "inout" or "harmonic" (centrality), default: "inout"} 
\item{legend}{should be colors and legend labels returned instead of colors, default: FALSE} 
#  } 
\examples{ 
par(mfrow=c(1,2),mai=rep(0.1,4)) 
A=mgraph(type="random",nodes=6,edges=8) 
cols=mgraph_nodeColors(A) 
mgraph_degree(A,mode="in") 
mgraph_degree(A,mode="out") 
cols 
plot(A, layout="star") 
plot(A, layout="star",vertex.color=cols) 
cols=mgraph_nodeColors(A,type="harmonic") 
plot(A, layout="star",vertex.color=cols) 
lg=mgraph_nodeColors(A,type="harmonic",legend=TRUE) 
legend("bottom",fill=lg$col, legend=lg$legend) 
} 
 
