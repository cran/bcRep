\name{sequences.junctionFrame}
\alias{sequences.junctionFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of junction frames of sequences
}
\description{
This function gives information about the proportion of in-frame and out-of-frame sequences.
}
\usage{
sequences.junctionFrame(data = NULL, relativeValues=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     Vector containing junction frame information, e.g. taken from IMGT/HighV-QUEST output 1_Summary(...).txt
}
  \item{relativeValues}{
     Shall relative or absolute values be returned? (default: TRUE)
}
 \item{...}{

}
}
\details{
Junction: coding region of the V-J or V-D-J junction from 2nd CYS (position 104) to J-PHE or J-TRP (position 118) in an IG or TR V-DOMAIN (\url{http://www.imgt.org/ligmdb/label})
}
\value{
Output is a data frame including proportions of in-frame and out-of-frame sequences, as well as proportion of sequences without any information about junction frames.
}
\references{
IMGT/LIGM-DB labels: \url{http://www.imgt.org/ligmdb/label}

IMGT/HighV-QUEST definition of functionality (and junction frame): \url{http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#func}
}
\author{
Julia Bischof
}

\seealso{
     \code{\link{sequences.functionality}}
}
\examples{
data(summarytab)
junctionfr<-sequences.junctionFrame(data = summarytab$JUNCTION_frame)
barplot(as.numeric(junctionfr),xlab="",xlim=c(0,2),legend.text = colnames(junctionfr),
     col=c("orange","darkblue","gray"),xaxt="n",main="Junction frame usage",ylab="proportion")
}


