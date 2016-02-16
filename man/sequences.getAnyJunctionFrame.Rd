\name{sequences.getAnyJunctionFrame}
\alias{sequences.getAnyJunctionFrame}
\alias{sequences.getInFrames}
\alias{sequences.getOutOfFrames}
\title{
Filter for in-frame or out-of-frame sequences
}
\description{
Filter IMGT/HighV-QUEST output files for in-frame or out-of-frame sequences. 

\code{sequences.getInFrames} filters only in-frame sequences

\code{sequences.getOutOfFrames} filters only out-of-frame sequences

\code{sequences.getAnyJunctionFrame} filters all sequences with any junction frame information (in-frame or out-of-frame)
}
\usage{
sequences.getAnyJunctionFrame(data = NULL)

sequences.getInFrames(data = NULL)

sequences.getOutOfFrames(data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     IMGT/HighV-QUEST output file with junction frame information, e.g. 1_Summary().txt, ...
}
}

\value{
Output is a data frame with the same columns like the input, but filtered for only in-frame or out-of-frame sequences or sequences having any junction frame information.
}
\references{
IMGT/LIGM-DB labels: \url{http://www.imgt.org/ligmdb/label}

IMGT/HighV-QUEST definition of functionality (and junction frame): \url{http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#func}

}
\author{
Julia Bischof
}

\seealso{
     \code{\link{sequences.getAnyJunctionFrame}}, \code{\link{sequences.getInFrames}}, \code{\link{sequences.getOutOfFrames}}, 
     
     \code{\link{sequences.junctionFrame}}, \code{\link{sequences.getAnyFunctionality}}, \code{\link{sequences.getProductives}}, 
     
     \code{\link{sequences.getUnproductives}}
}
\examples{
data(summarytab)
InFrameSequences<-sequences.getInFrames(summarytab)
OutOfFrameSequences<-sequences.getOutOfFrames(summarytab)
AnyJunctionFrSequences<-sequences.getAnyJunctionFrame(summarytab)
}

