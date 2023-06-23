\name{mgraph_accuracy} 
\alias{mgraph_accuracy} 
\title{quality measures for a predicted graph} 
\description{ 
   The function `mgraph_accuracy` measures the prediction 
   quality of a predicted graph in comparison to a known true one. 
   The graph comparisons are done on the basis of undirected graphs only. 
   Directed graphs are converted to undirected internally. 
 
   The most used accuracy measure to balance the different measures is the F-score where: 
 
\itemize{ 
   \item the F0.5 measure (beta=0.5) gives more weight on precision, less weight on sensitivity/recall. 
   \item the F1 measure (beta=1.0) balances the weight on precision and sensitivity/recall. 
   \item the F2 measure (beta=2.0) gives less weight on precision, more weight on sensitivity/recall 
  } 
} 
\usage{mgraph_accuracy(g.true,g.pred)} 
\arguments{ 
\item{g.true}{snha, mgraph or or adjacency matrix of a the true graph} 
\item{g.pred}{snha, mgraph or or adjacency matrix of a the predicted graph} 
} 
\value{returns list object with various accuracy measures, such as Sens(itivity), Spec(ificity), Prec(ision), Acc(uracy), Non-information-rate (NIR), balanced classification rate (BCR), F1 measure and Mathhews correlation coefficient (MCC)} 
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
 
