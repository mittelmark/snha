\name{mgraph_lms} 
\alias{mgraph_lms} 
\title{linear model backward variable selection} 
\description{ 
   This function is a simple version of linear model building using backward 
   variable selection. At each step the variable with the lowest decrease value 
   in the adjusted R-square value is removed from the model and the r-square value is stored in the return matrix. 
} 
\usage{mgraph_lms(x)} 
\arguments{ 
\item{x}{a data matrix or data frame} 
} 
\value{returns matrix with the adjusted r-square values attributed by the variables} 
\examples{ 
data(swiss) 
round(mgraph_lms(swiss)*100,2) 
W = mgraph(type='werner') 
data=t(snha_graph2data(W)) 
round(mgraph_lms(data)*100,2) 
} 
 
