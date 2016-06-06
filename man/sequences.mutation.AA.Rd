\name{sequences.mutation.AA}
\alias{sequences.mutation.AA}
\alias{plotSequencesMutationAA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect amino acid mutations in V-, FR- or CDR-regions
}
\description{
This functions detects amino acid mutations in V-, FR- or CDR-regions and returns a data frame containing proportions of mutations from amino acid x to amino acid y.
}
\usage{
sequences.mutation.AA(mutationtab = NULL, 
     sequence = c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"))
plotSequencesMutationAA(mutationAAtab = NULL, 
     showChange = c("no","hydropathy", "volume", "chemical"), title = NULL, PDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationtab}{
     IMGT output table 7_V-REGION-mutation-and-AA-change-table(...).txt
}
  \item{sequence}{
     Sequence region to be analyzed (V, FR1, FR2, FR3, CDR1, CDR2)
}
  \item{mutationAAtab}{
     Output table of sequences.mutation.AA()
}
  \item{showChange}{
     Shall hydropathy, volume, chemical or no changes be returned?
}
  \item{title}{
     Title of the plot
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
\code{sequences.mutation.AA()} returns a 20 x 20 data frame containing proportions of amino acid mutations in a given region. Original amino acids are in rows, mutated ones in columns.

\code{plotSequencesMutationAA()} returns a kind of heatmap containing proportion of amino acid changes and hydropathy, chemical or volume changes. Original amino acids are in rows, mutated ones in columns. The \code{PDF} character string should be only the project name (without ".pdf"). A figure called "PDF"_AA-mutation.pdf will be saved to the working directory.
}

\references{
IMGT/HighV-QUEST V-REGION mutation and AA change table: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#mut-table}

IMGT Index "Mutation": \url{http://www.imgt.org/IMGTindex/mutation.php}
}
\author{
Julia Bischof
}

\seealso{
\code{\link{sequences.mutation.AA}}, \code{\link{plotSequencesMutationAA}}, 

\code{\link{sequences.mutation}}, \code{\link{sequences.mutation.base}}, 

\code{\link{plotSequencesMutationBase}}
}
\examples{
data(mutationtab)
V.mut.AA<-sequences.mutation.AA(mutationtab = mutationtab, sequence = "V")
\dontrun{
plotSequencesMutationAA(mutationAAtab = V.mut.AA, showChange = "hydropathy")
}

}
