\name{plot.snha} 
 
\alias{plot.snha} 
\title{display network or correlation matrices of snha graphs} 
\description{The function `plot.snha` provides a simple display of network  
  graphs correlation matrices using  filled circles (vertices) to represent 
  variables and edges which connect the vertices with high absolute.  
  correlation values. Positive correlations are shown in black, negative 
  correlations are shown in red. For more information see the details 
  section. 
} 
\details{This is a plot function to display networks or correlation matrices 
         of 'snha' graph objects.  
         In case of bootstrapping the graph by using the `snha` function with the `prob=TRUE`  
         option lines in style full, broken and dotted lines are drawn if 
         they are found in more than 75, 50 or 25 percent of all re-samplings. 
         You can change these limits by using the `threshold` argument.} 
\usage{ 
\method{plot}{snha}( 
 x, 
 type = "network", 
 layout = "circle", 
 vertex.color = "salmon", 
 cex = 1, 
 vertex.size = 5, 
 edge.width = 2, 
 edge.color = c("grey70", "red"), 
 edge.text = NULL, 
 edge.cex = 0.8, 
 edge.pch = 0, 
 noise = FALSE, 
 hilight.chain = NULL, 
 chain.color = c("black", "red"), 
 star.center = NULL, 
 plot.labels = TRUE, 
 lty = 1, 
 threshold = c(0.25, 0.5, 0.75), 
 interactive = FALSE, 
 ... 
) 
} 
\arguments{ 
\item{x}{snha graph object usually created with the 'snha' function or an 
         adjacency matrix} 
\item{type}{character string specifying the plot type either 'network' or 
'         cor', default: 'network'} 
\item{layout}{graph layout for plotting one of 'circle', 'sam', 'samd', 
          'grid', 'mds', 'mdsd', 'star', default: 'circle'} 
\item{vertex.color}{default color for the vertices, either a single value, 
           all vertices have hen this color or a vector of values, 
           for different colors for the nodes, default: 'salmon'} 
\item{cex}{size of the vertex labels which are plotted on the vertices, default: 1} 
\item{vertex.size}{number how large the vertices should be plotted, default: 5} 
\item{edge.width}{number on how strong the edges should be plotted, if 
          edge.width=0, then the number is based on the correlation values,  
          default: 2} 
\item{edge.color}{color to be plotted for edges. Usually vector of length two. 
          First color for positive correlations, second color for negative 
          correlations. Default: c('grey','red')} 
\item{edge.text}{optional matrix to give edge labels, default: NULL} 
\item{edge.cex}{character expansion for edge labels, default: 0.8} 
\item{edge.pch}{plotting character which should be placed below the edge.text, 
           default: 0} 
\item{noise}{should be noise added to the layout. Sometimes useful if nodes 
           are too close. Default: FALSE} 
\item{hilight.chain}{which chain should be highlighted,  
           default: NULL (no chain highlight)} 
\item{chain.color}{which color for chain edges, default: black} 
\item{star.center}{the centered node if layout is 'start', must be 
           a character string for the node name, default: NULL} 
\item{plot.labels}{should node labels plotted, default: TRUE} 
\item{lty}{line type for standard edges in the graph, default: 1} 
\item{threshold}{cutoff values for bootstrap probabilities for drawing edges 
         as dotted, broken lines and solid lines, default: c(0.25,0.5,0.75)} 
\item{interactive}{switch into interactive mode where you can click in the 
         graph and move nodes with two clicks, first selecting the node, 
         second click gives the new coordinates for the node, default: FALSE} 
\item{\dots}{currently not used} 
} 
\value{returns the layout of the plotted network or NULL if type is `corrplot` 
        (invisible)} 
\examples{  
data(swiss) 
sw.g=snha(swiss,method='spearman') 
sw.g$theta 
round(sw.g$sigma,2) 
plot(sw.g,type='network',layout='circle') 
plot(sw.g,type='network',layout='sam') 
plot(sw.g,type='corplot') 
# adding correlation values 
plot(sw.g,edge.text=round(sw.g$sigma,2),edge.cex=1.2,edge.pch=15) 
sw.g=snha(swiss,method='spearman',prob=TRUE) 
sw.g$theta 
sw.g$probabilities 
plot(sw.g,type='network',layout='sam') 
sw.g$chains 
# plot chains for a node 
plot(sw.g,layout="sam",lty=2,hilight.chain="Infant.Mortality", 
 edge.width=3,edge.color=c("black","red")) 
# an example for an adjacency matrix 
M=matrix(rbinom(100,1, 0.2),nrow=10,ncol=10) 
diag(M)=0 
colnames(M)=rownames(M)=LETTERS[1:10] 
plot.snha(M) 
} 
 
