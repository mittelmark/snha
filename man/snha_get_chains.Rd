\name{snha_get_chains} 
\alias{snha_get_chains} 
\title{Return the chains of an snha graph as data frame} 
\usage{snha_get_chains(graph)} 
\description{This is a utility function to return the chains which 
constructs the graph as a matrix.} 
\arguments{ 
\item{graph}{a snha graph object} 
} 
\value{matrix with one chain per row, shorter chains are filled up with empty strings} 
\examples{ 
data(swiss) 
sw.g=snha(swiss) 
snha_get_chains(sw.g) 
} 
 
