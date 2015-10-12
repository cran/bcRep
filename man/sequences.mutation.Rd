\name{sequences.mutation}
\alias{sequences.mutation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Basic statistics on mutations of sequences
}
\description{
This function summarizes the number of mutations of sequences. It gives also information about the number of silent and replacement mutations, R/S ratio, as well as mutation numbers, depending on functionality or junction frame distributions.
}
\usage{
sequences.mutation(mutationtab = NULL, summarytab = NULL, 
     sequence = c("V", "FR1", "FR2", "FR3", "CDR1", "CDR2"), functionality = FALSE, 
     junctionFr = FALSE, rsRatio=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summarytab}{
     IMGT/HighV-QUEST output 1_Summarytab(...).txt
}
  \item{mutationtab}{
     IMGT/HighV-QUEST output 7_V-REGION-mutation-and-AA-change-table(...).txt
}
  \item{sequence}{
     One of V, FR1, FR2, FR3, CDR1, CDR2
}
  \item{functionality}{
     TRUE: mutation vs. functionality will be returned (default: FALSE)
}
  \item{junctionFr}{
     TRUE: mutation vs. junction frame usage will be returned, summarytab required (default: FALSE)
}
  \item{rsRatio}{
     TRUE: R/S ratio will be returned (default: FALSE)
}
 \item{...}{

}
}
\details{
IMGT/HighV-QUEST output 7_V-REGION-mutation-and-AA-change-table(...).txt (\code{mutationtab}) is required as input. 1_Summarytab(...).txt (\code{summarytab}) is optional; if specified, junction frame information and in special for V sequences "V-REGION identity [nt]" can be returned. Mutations of V region, as well as FR1, 2, 3 and CDR1, 2 can be analyzed.
\code{rsRatio=T} returns the ratio of replacement and silent mutations per sequence. Sequences, where there is either no silent or no replacement mutation, will have a ratio of 0.
}
\value{
Output is a list containing
  \item{Number_of_mutations}{data frame with number of total mutations, replacement and silent mutations (optional: V sequences V-REGION identity [nt], R/S ratio)}
  \item{Functionality}{Proportions of mutations and no mutations in productive and unproductive sequences (optional)}
  \item{Junction_frame}{Proportions of mutations and no mutations in in-frame and out-of-frame sequences (optional)}
    \item{RS_ratio}{Ratio of replacement and silent mutations (optional)}
}
\references{
IMGT/HighV-QUEST V-REGION mutation and AA change table: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#mut-table}

IMGT Index "Mutation": \url{http://www.imgt.org/IMGTindex/mutation.html}
}
\author{
Julia Bischof
}

\examples{
data(mutationtab)
data(summarytab)
V.mutation<-sequences.mutation(mutationtab = mutationtab, summarytab = summarytab, 
     sequence = "V", junctionFr = TRUE, rsRatio=TRUE)
CDR1.mutation<-sequences.mutation(mutationtab = mutationtab, sequence = "CDR1",
     functionality=TRUE)
par(mar=c(18,5,5,3))
barplot(as.numeric(CDR1.mutation$Functionality[,1]),
     names=rownames(CDR1.mutation$Functionality),
     ylab="proportion",main="Mutation vs. Functionality",las=3)
}



