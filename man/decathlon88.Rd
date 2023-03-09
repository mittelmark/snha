\docType{data}
\name{decathlon88}
\alias{decathlon88}
\title{Men Decathlon data from the 1988 Olympics}
\usage{decathlon88}
\format{
A data frame with 33 rows and 10 columns:
\describe{
  \item{disc}{discus results in m}
  \item{high}{high jump results in m}
  \item{jave}{javelin through results in m}
  \item{long}{long jump results in m}
  \item{pole}{pole vault results in m}
  \item{shot}{shot put results in m}
  \item{X100}{running speed over 100m in km/h}
  \item{X110}{running speed over 110m hurdles in km/h}
  \item{X1500}{running speed over 1500m in km/h}
  \item{X400}{running speed over 400m in km/h}
}
}
\source{
<https://en.wikipedia.org/wiki/Athletics_at_the_1988_Summer_Olympics_-_Men's_decathlon>
}
\description{
A subset of data from the Decathlon from the 1988 Olympic games. 
Included are all athletes which finished with more than 7000 points.
}
\examples{
data(decathlon88)
head(decathlon88)
A=snha(decathlon88,method="spearman",alpha=0.1)
cols=rep("salmon",10)
cols[names(A$data) \%in\% c("jave","shot","disc","pole")]="skyblue"
plot(A,layout="sam",vertex.color=cols,vertex.size=8,cex=1.2,edge.width=5)
snha_rsquare(A)

}
\keyword{datasets}
