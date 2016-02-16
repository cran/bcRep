\name{clones.filterJunctionFrame}
\alias{clones.filterJunctionFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter for in-frame/out-of-frame clones
}
\description{
Clones containing only in-frame or only out-of-frame sequences can be filtered.
}
\usage{
clones.filterJunctionFrame(clones.tab = NULL, 
     filter = c("in-frame", "out-of-frame"))
}
\arguments{
  \item{clones.tab}{
     A data frame containing clones and their characteristics
}
  \item{filter}{
     Filter feature
}
}
\value{
Output is a data frame, that contains clones with only in-frame or out-of-frame sequences.

}
\references{
IMGT/LIGM-DB labels: \url{http://www.imgt.org/ligmdb/label}

IMGT/HighV-QUEST definition of functionality (and junction frame): \url{http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#func}
}
\author{
Julia Bischof
}
\seealso{
\code{\link{clones.filterFunctionality}}, \code{\link{clones.filterSize}}
}

\examples{
data(clones.ind)

inFrameClones<-clones.filterJunctionFrame(clones.tab = clones.ind, 
     filter = "in-frame")
}
