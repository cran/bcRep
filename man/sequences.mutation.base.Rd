\name{sequences.mutation.base}
\alias{sequences.mutation.base}
\alias{plotSequencesMutationBase}

\title{
Statistics about silent mutations
}
\description{
This function calculates a) proportions of silent mutations from nucleotide A to B and b) proportions of A, C, G and T/U from positions -3 to +3 around a mutation.
}
\usage{
sequences.mutation.base(mutationtab = NULL, summarytab = NULL, 
     analyseEnvironment = FALSE, analyseMutation = TRUE,
     sequence = c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"), nrCores = 1)


plotSequencesMutationBase(mutationBaseTab = NULL, plotEnvironment = FALSE, 
     plotMutation = TRUE, colHeatmap = c("white","darkblue"), title = NULL, PDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationtab}{
     IMGT output table 7_V-REGION-mutation-and-AA-change-table(...).txt
}
  \item{summarytab}{
     IMGT output table 1_Summary(...).txt
}
 \item{analyseEnvironment}{
     Shall proportions of A, C, G and T/U from positions -3 to +3 around a mutation be analysed? (default: FALSE)
}
 \item{analyseMutation}{
     Shall proportions of mutations from nucleotide A to B be analysed? (default: TRUE)
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
 \item{plotEnvironment}{
     Shall proportions of A, C, G and T/U from positions -3 to +3 around a mutation be plotted? (default: FALSE)
}
 \item{plotMutation}{
     Shall proportions of mutations from nucleotide A to B be plotted? (default: TRUE)
}
 \item{colHeatmap}{
     Colours for the heatmap (default: c("white","darkblue"))
}
\item{title}{
     Title of the plot
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
This function calculates a) proportions of silent mutations from nucleotide A to B and b) proportions of A, C, G and T/U from positions -3 to +3 around a mutation. For case a) a matrix containing proportions of mutations is returned. In case b) a data frame containing positions in columns and base changes in rows is returned. Position 0 is the position where mutation occurs. Position -3 to -1 and +1 to +3 show proportions of bases next to the mutated position.


Function \code{plotSequencesMutationBase()} returns a heatmap or barplots for every base. The \code{PDF} character string should be only the project name (without ".pdf"). A figure called "PDF"_Base-mutation_mutated-position.pdf" and/or "PDF"_Base-mutation_environment.pdf" will be saved to the working directory.
}
\references{
IMGT/HighV-QUEST V-REGION mutation and AA change table: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#mut-table}

IMGT Index "Mutation": \url{http://www.imgt.org/IMGTindex/mutation.php}
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
plotSequencesMutationBase(mutationBaseTab = V.base.mut, plotMutation = T)
}

}

