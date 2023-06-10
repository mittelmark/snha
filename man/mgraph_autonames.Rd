\name{mgraph_autonames} 
\alias{mgraph_autonames} 
\title{create names for nodes and other data structures} 
\description{ 
    This function aids in creating standard node labels for graphs. 
} 
\usage{mgraph_autonames( 
  n, 
  prefix=NULL) 
} 
\arguments{ 
\item{n}{how many labels} 
\item{prefix}{the node prefix, per defaultif not given LETTERS with numbers will be used, default: NULL} 
#  } 
\examples{ 
mgraph_autonames(12) 
mgraph_autonames(12,LETTERS[1:4]) 
mgraph_autonames(12,"R") 
} 
 
