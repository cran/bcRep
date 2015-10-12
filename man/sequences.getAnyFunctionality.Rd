\name{sequences.getAnyFunctionality}
\alias{sequences.getAnyFunctionality}
\alias{sequences.getProductives}
\alias{sequences.getUnproductives}
\title{
Filter for productive or unproductive sequences
}
\description{
Filter IMGT/HighV-QUEST output files for productive or unproductive sequences. 

\code{sequences.getProductives} filters only productive sequences

\code{sequences.getUnproductives} filters only unproductive sequences

\code{sequences.getAnyFunctionality} filters all sequences with any functionality information (productive or unproductive)
}
\usage{
sequences.getAnyFunctionality(data = NULL)

sequences.getProductives(data = NULL)

sequences.getUnproductives(data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     IMGT/HighV-QUEST output file with functionality information, f.e. 1_Summary().txt, 3_Nt-sequences.txt, ...
}
}

\value{
Output is a data frame with the same columns like the input, but filtered for only productive, unproductive or sequences having any functionality information.
}
\references{
IMGT/LIGM-DB labels: \url{http://www.imgt.org/ligmdb/label}

IMGT/HighV-QUEST definition of functionality: \url{http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#func}
}
\author{
Julia Bischof
}

\seealso{
     \code{\link{sequences.getAnyFunctionality}}, \code{\link{sequences.getProductives}}, \code{\link{sequences.getUnproductives}}, 
     
     \code{\link{sequences.functionality}}, \code{\link{sequences.getAnyJunctionFrame}}, \code{\link{sequences.getInFrames}}, 
     
     \code{\link{sequences.getOutOfFrames}}
}
\examples{
data(summarytab)
ProductiveSequences<-sequences.getProductives(summarytab)
UnroductiveSequences<-sequences.getUnproductives(summarytab)
AnyFunctSequences<-sequences.getAnyFunctionality(summarytab)
}

