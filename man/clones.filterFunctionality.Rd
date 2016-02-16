\name{clones.filterFunctionality}
\alias{clones.filterFunctionality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter for productive/unproductive clones
}
\description{
Clones containing only productive or only unproductive sequences are filtered.
}
\usage{
clones.filterFunctionality(clones.tab = NULL, 
     filter = c("productive", "unproductive"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clones.tab}{
     A data frame containing clones and their characteristics
}
  \item{filter}{
     Filter feature
}
}

\value{
Output is a data frame, that contains clones with only productive or unproductive sequences. 
}
\references{
IMGT/LIGM-DB labels: \url{http://www.imgt.org/ligmdb/label}

IMGT/HighV-QUEST definition of functionality: \url{http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html#func}

}
\author{
Julia Bischof
}
\seealso{
\code{\link{clones.filterJunctionFrame}}, \code{\link{clones.filterSize}}
}
\examples{
data(clones.ind)

productiveClones<-clones.filterFunctionality(clones.tab = clones.ind, 
     filter = "productive")
}

