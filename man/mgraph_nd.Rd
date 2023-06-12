\name{mgraph_nd} 
\alias{mgraph_nd} 
\title{network deconvolved data matrix using the algorithm of Feizi et. al. 2013} 
\description{ 
   This function is a R implementation for "A General Method to Distinguish Direct Dependencies over Networks" by Feizi et. al (2013). For the Matlab code look here: 
  \url{http://compbio.mit.edu/nd/code/ND.m} 
} 
\usage{mgraph_nd(x,beta=0.99,alpha=1,control=FALSE)} 
\arguments{ 
\item{x}{symmetric relevance matrix, where high similarity between nodes is expressed with high values such as in an correlation matrix} 
\item{beta}{scaling paramater, the largest absolute eigenvalue is mapped to beta, values should be between 0 and 1, default: 0.99} 
\item{alpha}{fraction of edges of the observed dependency matrix to be kept, should be between 0 (none) and 1 (all), default: 1} 
\item{control}{if FALSE only direct weights are displayed, if TRUE also non-observed interactions are displayed, default: FALSE} 
} 
\value{returns matrix with deconvoluted relevance matrix} 
\examples{ 
W=mgraph(type="werner") 
W 
data=snha_graph2data(W) 
C=cor(t(data)) 
round(C,2) 
round(mgraph_nd(C,beta=0.9,alpha=0.3),2) 
} 
\references{ 
\itemize{ 
\item Feizi, S., Marbach, D., Medard, M., & Kellis, M. (2013). Network deconvolution as a general method to distinguish direct dependencies in networks. 
    Nature biotechnology, 31(8), 726-733. DOI 10.1038/nbt.2635 
} 
} 
\author{ 
\itemize{ 
\item @2013 KELLIS-LAB, Soheil Feizi, Matlab code 
\item @2021 Detlef Groth, University of Potsdam, R code 
} 
} 
 
