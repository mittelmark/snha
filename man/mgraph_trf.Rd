\name{mgraph_trf} 
\alias{mgraph_trf} 
\title{checking possible chains with 3 nodes for triad structures} 
\description{ 
   Experimental: 
 
   This function checks a given snha graph for possible triad structures 
   which were not recognized by the original algorithm. 
   The method checks all nodes within path length 2 for a possible edge 
   by comparing the ratio of the correlation coefficients. 
   The function snha uses defaults of 0.02 for del(eting) and 0.05 for add(ing) values. 
} 
\usage{mgraph_trf(x,frac=1.2)} 
\arguments{ 
\item{x}{snha graph object} 
\item{frac}{amount which the absolute correlation exceeds in compariuson of the product of the other two correlations in the possible triad, value should be larger than 1, default: 1.2} 
} 
\value{returns snha graph object with added triads to chains and to the adjacency matrix} 
\examples{ 
set.seed(124) 
W=mgraph(type="werner") 
data=t(snha_graph2data(W)) 
aa=snha(data) 
aa$theta 
unlist(mgraph_accuracy(W,aa$theta)) 
# remove false positive by linear model 
ab=mgraph_lmc(aa) 
aa$theta=ab 
aa=mgraph_trf(aa) 
unlist(mgraph_accuracy(W,aa$theta)) 
} 
 
