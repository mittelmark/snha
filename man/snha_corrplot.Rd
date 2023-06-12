\name{snha_corrplot} 
\alias{snha_corrplot} 
\title{visualize a matrix of pairwise correlations} 
\usage{snha_corrplot( 
   x, 
   text.lower=FALSE, 
   text.upper=FALSE, 
   pch.minus=19, 
   pch.plus=19, 
   xtext=NULL, 
   cex=1.0, 
   ...) 
} 
\description{This function returns xy coordinates for a given input 
      adjacency matrix or snha graph. It is useful if you like to plot 
      the same set of nodes with different edge connections  
      in the same layout.} 
\arguments{ 
 \item{x}{matrix with pairwise correlations} 
 \item{text.lower}{should in the lower diagonal the correlation coefficient be shown, default: TRUE} 
 \item{text.upper}{should in the upper diagonal the correlation coefficient be shown, default: FALSE} 
 \item{pch.minus}{the plotting symbol for negative correlations, default: 19} 
 \item{pch.plus}{the plotting symbol for positive correlations, default: 19} 
 \item{xtext}{labels which should be placed at the bottom of the plot,  default: NULL} 
 \item{cex}{character expansion for text and correlation symbols, default: 1} 
 \item{\ldots}{arguments delegated to the plot function} 
} 
\examples{ 
data(swiss) 
sw=swiss 
colnames(sw)=abbreviate(colnames(swiss),6) 
options(warn=-1) # avoid spearman warnings 
snha_corrplot(cor(sw,method="spearman"),cex.sym=8,text.lower=TRUE) 
options(warn=0) 
} 
 
