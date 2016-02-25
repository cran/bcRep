\name{sequences.mutation.base}
\alias{sequences.mutation.base}
\alias{plotSequencesMutationBase}

\title{
Base proportion nearby mutations
}
\description{
This function calculates proportion of A, C, G and T from positions -3 to +3 around a mutation.
}
\usage{
sequences.mutation.base(mutationtab = NULL, summarytab = NULL, 
     sequence = c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"), nrCores = 1)


plotSequencesMutationBase(mutationBaseTab = NULL, title = NULL, PDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationtab}{
     IMGT output table 7_V-REGION-mutation-and-AA-change-table(...).txt
}
  \item{summarytab}{
     IMGT output table 1_Summary(...).txt
}
  \item{sequence}{
     Sequence region to be analyzed (V, FR1, FR2, FR3, CDR1, CDR2)
}
  \item{nrCores}{
     Number of cores to be used for parallel processing
}
  \item{mutationBaseTab}{
     Output of sequences.mutation.base()
}
  \item{title}{
     Title of the plot
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
This function calculates proportion of A, C, G and T from positions -3 to +3 around a mutation. Position 0 is the position where mutation occurs. Position -3 to -1 and +1 to +3 show proportions of bases next to the mutated position.


Function \code{plotSequencesMutationBase()} returns barplots for every base. The \code{PDF} character string should be only the project name (without ".pdf"). A figure called "PDF"_Base-mutation.pdf will be saved to the working directory.
}
\value{
\code{sequences.mutation.base()} returns a data frame containing positions in columns and base changes in rows.
}
\references{
IMGT/HighV-QUEST V-REGION mutation and AA change table: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#mut-table}

IMGT Index "Mutation": \url{http://www.imgt.org/IMGTindex/mutation.html}
}
\author{
Julia Bischof
}

\seealso{
\code{\link{sequences.mutation.base}}, \code{\link{plotSequencesMutationBase}}, 

\code{\link{sequences.mutation.AA}}, \code{\link{plotSequencesMutationAA}}, 

\code{\link{sequences.mutation}}, 
}
\examples{
data(mutationtab)
data(summarytab)
V.base.mut<-sequences.mutation.base(mutationtab = mutationtab, summarytab = summarytab, 
     sequence = "V", nrCores = 1)

\dontrun{
plotSequencesMutationBase(mutationBaseTab = V.base.mut)
}

}

