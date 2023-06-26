\name{snha} 
\alias{snha} 
\title{Initialize a snha object with data.} 
\description{ 
    The main entry function to initialize a snha object with data where 
    variables are in columns and items are in rows. If you like to use the network 
    deconvolution method of Feizi et al. instead of the data submit the deconvoluted matrix. 
    As methods either correlation methods such as pearson, spearman or kendall or mutual information (mi) can be used. 
   The latter method as  well  covers non-linear linear relationships. 
} 
\usage{snha( 
  data, 
  alpha=0.05, 
  method='pearson', 
  threshold=0.01, 
  check.singles=FALSE, 
  prob=FALSE, 
  prob.threshold=0.2, 
  prob.n=25, 
  nd=FALSE, 
  lmc=FALSE, 
  lma=FALSE) 
} 
\arguments{ 
\item{data}{a data frame where network nodes are the row names and data 
            variables are in the columns.} 
\item{alpha}{confidence threshold for p-value edge cutting after all chains 
             were generated, default: 0.05.} 
\item{method}{method to calculate correlation/association values, can be 
              'pearson', 'spearman' or 'kendall', 'mi', default: 'pearson'.} 
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
\item{nd}{perfom network deconvolution due to the method of Feizi et. al. 2013, default: FALSE} 
\item{lmc}{perfom a linear model check, edges which does not improve the node model by adjacent edges more than 0.02 are removed, default: FALSE} 
\item{lma}{perfom a linear model check, edges which does improve the node model by adjacent edges with more than 0.02 R-square are added, works only if lmc is TRUE,default: FALSE} 
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
sw.m=snha(swiss,method="mi") 
sw.m$theta 
par(mfrow=c(1,2)) 
plot(sw.g) 
plot(sw.m) 
} 
\seealso{ \link[snha:plot.snha]{plot.snha}} 
 
