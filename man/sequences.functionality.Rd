\name{sequences.functionality}
\alias{sequences.functionality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of functionality of sequences
}
\description{
This function gives information about the proportion of productive and unproductive sequences.
}
\usage{
sequences.functionality(data = NULL, relativeValues=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     Vector containing functionality information, f.e. taken from IMGT/HighV-QUEST output, f.e. 1_Summary(...).txt, 3_Nt-sequences(...).txt, 5_AA-sequences(...).txt
}
  \item{relativeValues}{
     Shall relative or absolute values be returned? (default: TRUE)
}
 \item{...}{

}
}

\value{
Output is a data frame including proportions of productive and unproductive sequences, as well as proportion of sequences with unknown functionality.
}
\references{
IMGT/LIGM-DB labels: \url{http://www.imgt.org/ligmdb/label}

IMGT/HighV-QUEST definition of functionality: \url{http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#func}
}
\author{
Julia Bischof
}
\seealso{
     \code{\link{sequences.junctionFrame}}
}
\examples{
data(summarytab)
funct<-sequences.functionality(data = summarytab$Functionality)
barplot(as.numeric(funct),xlab="",xlim=c(0,2),legend.text = colnames(funct),
     col=c("orange","darkblue","gray"),xaxt="n",main="Functionality",ylab="proportion")
}


