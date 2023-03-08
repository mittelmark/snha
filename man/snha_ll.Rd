\name{snha_ll} 
\alias{snha_ll} 
\title{log-likelihood for the given snha graph and the given chain} 
\description{This function returns the log-likelihood for the given snha 
graph and the given chain.} 
\arguments{ 
\item{g}{a snha graph object} 
\item{chain}{a chain object of a snha graph, if not given a data frame with 
the values is returned for all chains, default: NULL} 
} 
\value{list with the following components: 'll.total', 'll.chain', 
'll.rest', 'll.block', data frame 'df' with the columns 'chisq', 'p.value', 
'block.df', 'block.ch', 'block.p.value'. If chain is not given an overall 
summary is made for all chains an returned as data frame.} 
\examples{ 
data(swiss) 
sw.g=snha(swiss) 
snha_ll(sw.g,sw.g$chain$Catholic) 
head(snha_ll(sw.g)) 
} 
 
