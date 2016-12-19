\name{clones.IDlist}
\alias{clones.IDlist}
\docType{data}
\title{
Match sequence ID's and clone ID's
}
\description{
Get information about which sequences belong to the same clone (see details). This analysis can only be performed, if sequence ID's are returned in table containing clones (\code{clones(..., dispSeqID = T)}).
}
\usage{clones.IDlist(clones.seqID = NULL, summarytab.seqID = NULL)
}
\arguments{
  \item{clones.seqID}{
	A vector containing the column "Sequence_IDs" of the output of function clones()
  }
    \item{summarytab.seqID}{
	A vector containing the column "Sequence_ID" of 1_Summary(...).txt IMGT table
  }
}
\details{
This function returns information about each sequence given in den IMGT tables and the corresponding clone. This analysis can only be performed, if sequence ID's are returned in table containing clones (\code{clones(..., dispSeqID = T)}).

There are two columns: Sequence_ID and Clone_ID, where Sequence_ID is the same ID given in IMGT tables and Clone_ID is a ID representing all the given clones (clone_1, ..., clone_n, where n is the number of clones; "no_clone" represents the sequences, which belong to none of the clones [probably due to too small abundances]).

If \code{summarytab.seqID} is provided, the output table will be ordered like the sequences of the IMGT tables; including those sequences, which belong to no clone. If not, sequences are ordered for sequence ID of the \code{clones.seqID} vector.
}
\author{
Julia Bischof
}
\seealso{
 \code{\link{clones}}
}
\examples{
\dontrun{
data(summarytab)
data(aaseqtab)
clone.table<-clones(summarytab = summarytab, aaseqtab = aaseqtab, useJ = T, 
                         dispSeqID = T)

clone.ID<-clones.IDlist(clones.seqID = clone.table$Sequence_IDs, 
                        summarytab.seqID = summarytab$Sequence_ID)
}
}
