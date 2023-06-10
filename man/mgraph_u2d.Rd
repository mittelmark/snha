\name{mgraph_u2d} 
\alias{mgraph_u2d} 
\title{reate a directed graph out of an undirected graph} 
\description{ 
  This function creates a directed graph from an undirected one by the given input nodes. Input nodes can be chosen by names or a number for random selection of input nodes will be given. 
Input nodes will have at least shortest path distance to other input nodes of pathlength two. 
Selected input nodes will draw in each iteration outgoing edges to other nodes in the nth iteration neighborhood. The input 
nodes will alternatively select the next edges on the path to not visited nodes. All edges will be only visited onces. 
} 
\usage{mgraph_u2d( 
    g,input=2,negative=0.0,shuffle=FALSE 
)} 
\arguments{ 
\item{g}{mgraph object created with mgraph} 
\item{input}{number or names of input nodes in the graph, if number of input nodes is smaller than number of components, for each component one input node is automatically created} 
\item{negative}{proportion of inhibitive associations in the network value between 0 and 1 are acceptable, Default 0.0} 
\item{shuffle}{should just the edge directions beeing shuffled, if TRUE the graph will be very random without a real structure or chains of associations, default: FALSE} 
#  } 
\examples{ 
par(mfrow=c(2,2),mai=rep(0.1,4)) 
G=mgraph(type="angie",nodes=7,edges=9) 
mgraph_degree(G,mode='in') 
U=mgraph_d2u(G) 
H=mgraph_u2d(U,input=c("G","E")) 
identical(G,H) 
I=mgraph_u2d(U,shuffle=TRUE) 
identical(G,I) 
lay=snha_layout(G) 
plot(G,layout=lay,vertex.color=mgraph_nodeColors(G)) 
plot(U,layout=lay) 
plot(H,layout=lay,vertex.color=mgraph_nodeColors(H)) 
plot(I,layout=lay,vertex.color=mgraph_nodeColors(I)) 
} 
 
