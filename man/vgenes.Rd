\name{vgenes}
\alias{vgenes}
\docType{data}
\title{
VH gene usage data
}
\description{
VH gene usage data of 10 samples (rows) and 30 genes (columns).
}
\usage{data("vgenes")}
\format{
  The format is:
 num [1:10, 1:30] 0.0628 0.0529 0.0248 0.0512 0.0184 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:10] "Sample1" "Sample2" "Sample3" "Sample4" ...
  ..$ : chr [1:30] "IGHV4-34" "IGHV1-8" "IGHV1-3" "IGHV3-23" ...
}
\examples{
data(vgenes)
str(vgenes)
barplot(t(vgenes), col = rainbow(n = ncol(vgenes),start=0.2, end = 0.9), 
     xlim = c(0, nrow(vgenes)+5), ylab="proportion", 
     main="VH gene usage", las=3)
legend("right",col = rainbow(n = ncol(vgenes)), colnames(vgenes), pch=15)
}
\keyword{datasets}
