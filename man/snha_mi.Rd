\name{snha_mi} 
\alias{snha_mi} 
\title{mutual information for two vectors or a matrix} 
\usage{snha_mi(x,y=NULL,breaks=4,coc=FALSE,gnorm=FALSE)} 
\description{This function computes the mutual information between two numerical variables 
 or for all pairs of variables if a matrix is given. The normalization is actual value divided by the the geometric mean: I'(X,Y) = I(X,Y)/sqrt(I(X)*I(Y))  
} 
\arguments{ 
\item{x}{either a binned table, a matrix or data.frame or a numerical vector} 
\item{y}{numerical vector, required if x is a vector, otherwise ignored,default:NULL} 
\item{breaks}{number of breaks to create a binned table, default: 4} 
\item{coc}{Coefficient of constrains, Cxy=Ixy/Hy and Cyx=Ixy/Hx, this measure is asymetric,default:FALSE} 
\item{gnorm}{mutual information should be normalized by using the geometric mean Ixy=Ixy/sqrt(Ix*Iy),default:FALSE} 
} 
\value{mutual information value as scalar if input is table or two vectors or as matrix if input is matrix or data.frame} 
\examples{ 
rn1=rnorm(100,mean=10,sd=1); 
rn2=rn1+0.5*rnorm(100) 
cor(rn1,rn2) # high 
cor(rn1,sample(rn2)) #random  
snha_mi(rn1,rn2) # high  
snha_mi(rn1,sample(rn2)) #random 
snha_mi(rn1,rn2,breaks=4) 
snha_mi(rn1,rn2,breaks=7) 
data(swiss) 
round(snha_mi(swiss),2) 
round(snha_mi(swiss,gnorm=TRUE),2) 
as=snha(snha_mi(swiss,gnorm=TRUE)) 
plot(as) 
} 
