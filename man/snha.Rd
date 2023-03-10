\name{snha} 
\alias{snha} 
\title{Initialize a snha object with data.} 
\description{ 
    The main entry function to initialize a snha object with data where 
    variables are in columns and items are in rows 
} 
\usage{snha( 
  data, 
  alpha=0.05, 
  method='pearson', 
  threshold=0.01, 
  check.singles=FALSE, 
  prob=FALSE, 
  prob.threshold=0.2, 
  prob.n=25) 
} 
\arguments{ 
\item{data}{a data frame where network nodes are the row names and data 
            variables are in the columns.} 
\item{alpha}{confidence threshold for p-value edge cutting after all chains 
             were generated, default: 0.05.} 
\item{method}{method to calculate correlation/association values, can be 
              'pearson', 'spearman' or 'kendall', default: 'pearson'.} 
\item{threshold}{R-squared correlation coefficient threshold for which 
          r-square values should be used for chain generation, r=0.1 is r-square of 
          0.01, default: 0.01.} 
\item{check.singles}{should isolated nodes connected with sufficient high R^2 
          and significance, default: FALSE.} 
\item{prob}{should be probabilities computed for each edge using 
         bootstrapping. Only in this case the parameters starting with prob are used, 
         default: FALSE} 
\item{prob.threshold}{threshold to set an edge, a value of 0.5 means, that 
         the edge must be found in 50\% of all samplings, default: 0.2} 
\item{prob.n}{number of bootstrap samples to be taken, default: 25} 
} 
\keyword{network} 
\keyword{correlation} 
\value{A snha graph data object with the folling components: 
\describe{ 
\item{chains}{association chains building the graph} 
\item{data}{representing the original input data} 
\item{p.values}{matrix with p-values for the pairwise correlations} 
\item{probabilities}{in case of re-samplings, the proportion how often the chain was found} 
\item{sigma}{correlation matrix used for the algorithm} 
\item{theta}{adjacency  matrix found by the SNHA method} 
} 
} 
\examples{ 
data(swiss) 
sw.g=snha(swiss,method='spearman') 
# what objects are there? 
ls(sw.g) 
sw.g$theta 
round(sw.g$sigma,2) 
sw.g=snha(swiss,method='spearman',check.singles=TRUE,prob=TRUE) 
sw.g$theta 
sw.g$probabilities 
} 
\seealso{ \link[snha:plot.snha]{plot.snha}} 
 
