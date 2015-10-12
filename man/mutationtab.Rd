\name{mutationtab}
\alias{mutationtab}
\docType{data}
\title{
Summary of mutations cells 
}
\description{
Data frame represents an output from IMGT/HighV-QUEST, file 7_V-REGION-mutation-and-AA-change-table(...).txt}
\usage{data("mutationtab")}
\format{
  A data frame with 3000 observations on the following 11 variables.
  \describe{
    \item{\code{Sequence_number}}{a numeric vector}
    \item{\code{Sequence_ID}}{a character vector}
    \item{\code{Functionality}}{a character vector}
    \item{\code{V_GENE_and_allele}}{a character vector}
    \item{\code{V_REGION}}{a character vector}
    \item{\code{FR1_IMGT}}{a character vector}
    \item{\code{CDR1_IMGT}}{a character vector}
    \item{\code{FR2_IMGT}}{a character vector}
    \item{\code{CDR2_IMGT}}{a character vector}
    \item{\code{FR3_IMGT}}{a character vector}
    \item{\code{CDR3_IMGT}}{a character vector}
  }
}

\references{
Alamyar, E. et al., IMGT/HighV-QUEST: A High-Throughput System and Web Portal for the Analysis of Rearranged Nucleotide Sequences of Antigen Receptors, JOBIM 2010, Paper 63 (2010).

Brochet, X. et al., Nucl. Acids Res. 36, W503-508 (2008). PMID: 18503082

IMGT/HighV-QUEST 7_V-REGION-mutation-and-AA-change-table file description: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#Emuttable}
}
\examples{
data(mutationtab)
str(mutationtab)
}
\keyword{datasets}
