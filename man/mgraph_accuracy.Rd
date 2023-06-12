\name{mgraph_accuracy} 
\alias{mgraph_accuracy} 
\title{quality measures for a predicted graph} 
\description{ 
   The function `mgraph_accuracy` measures the prediction 
   quality of a predicted graph in comparison to a known true one. 
   The graph comparisons are done on the basis of undirected graphs only. 
   Directed graphs are converted to undirected internally. 
} 
\usage{mgraph_accuracy(g.true,g.pred)} 
\arguments{ 
\item{g.true}{snha, mgraph or or adjacency matrix of a the true graph} 
\item{g.pred}{snha, mgraph or or adjacency matrix of a the predicted graph} 
} 
\value{returns list object with various accuracy measures, such as Sens(itivity), Spec(ificity), Acc(uracy), balanced classification rate (BCR), F1 measure and Mathhews correlaition coefficient (MCC)} 
\examples{ 
ang=mgraph(type="angie",nodes=12,edges=18) 
data=snha_graph2data(ang) 
pred=snha(t(data))$theta 
round(unlist(mgraph_accuracy(ang,pred)),2) 
#  using vectors 
tvec=c(1,1,1,0,0,0) 
pvec=c(1,1,0,1,0,0) 
table(tvec,pvec) 
mgraph_accuracy(tvec,pvec) 
} 
 
