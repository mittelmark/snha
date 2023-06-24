\name{mgraph_lmc} 
\alias{mgraph_lmc} 
\title{using linear models to check the create snha graph} 
\description{ 
   This function checks a given snha graph against its own data using linear 
   models. Edges are removed if they do not add  more than 2 percent of the  
   adjusted  R-square value. The function might be only called directly if the user would like to 
   change the default values for del and add.  
   The function snha uses defaults of 0.02 for del(eting) and 0.05 for add(ing) values. 
} 
\usage{mgraph_lmc(x,del=0.02,add=NULL)} 
\arguments{ 
    \item{x}{snha graph object} 
    \item{del}{r-square threshold to delete edges, edges not adding more than this value to the linear model to the target node will be deleted, default: 0.02} 
    \item{add}{r-square threshold to add edges, if NULL no edges will be added, recommended value is 0.05, default: NULL} 
} 
\value{returns adjacency matrix where edges migh be removed if they are not adding explanation to the model} 
\examples{ 
set.seed(123) 
B=mgraph(type="barabasi",m=2,nodes=6) 
data=t(snha_graph2data(B)) 
aa=snha(data) 
aa$theta 
unlist(mgraph_accuracy(B,aa$theta)) 
ab=mgraph_lmc(aa) 
unlist(mgraph_accuracy(B,ab)) 
ab 
aa$theta-ab 
ac=mgraph_lmc(aa,add=0.02) 
unlist(mgraph_accuracy(B,ac)) 
} 
 
