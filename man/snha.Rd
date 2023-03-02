\name{snha} 
\alias{snha} 
\title{Initialize a snha object with data.} 
\description{ 
    The main entry function to initialize a snha object with data where 
    variables are in columns and items are in rows 
} 
\usage{ snha(x, method="pearson") } 
\arguments{ 
  \item{x}{ 
    data frame or matrix, variables in columns, items in rows. 
  } 
  \item{method}{ 
    The method to determine pairwise associations between variables such as 
    'pearson' or 'spearman', default: 'pearson' 
  } 
} 
\details{ 
    Some more details ... 
} 
\value{return an object of class snha with the following components: 
\begin{itemize}{ 
\item data - the input data 
\item sigma - the correlation matrix 
} 
} 
\examples{ 
    library(MASS) 
    data(birthwt) 
    head(birthwt) 
    as=snha(birtwth) 
    ls(as) 
} 
\seealso{  See also \code{\link{testprint}} } 
\keyword{ arith } 
 
