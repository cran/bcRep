\name{summarytab}
\alias{summarytab}
\docType{data}
\title{
Summary information of B cells of one individual
}
\description{
Data frame represents an output from IMGT/HighV-QUEST, file 1_Summary(...).txt
}
\usage{data("summarytab")}
\format{
  A data frame with 3000 observations on the following 29 variables.
  \describe{
    \item{\code{Sequence_number}}{a numeric vector}
    \item{\code{Sequence_ID}}{a character vector}
    \item{\code{Functionality}}{a character vector}
    \item{\code{V_GENE_and_allele}}{a character vector}
    \item{\code{V_REGION_score}}{a numeric vector}
    \item{\code{V_REGION_identity_\%}}{a numeric vector}
    \item{\code{V_REGION_identity_nt}}{a character vector}
    \item{\code{V_REGION_identity_\%_(with_ins/del_events)}}{a numeric vector}
    \item{\code{V_REGION_identity_nt_(with_ins/del_events)}}{a character vector}
    \item{\code{J_GENE_and_allele}}{a character vector}
    \item{\code{J_REGION_score}}{a numeric vector}
    \item{\code{J_REGION_identity_\%}}{a numeric vector}
    \item{\code{J_REGION_identity_nt}}{a character vector}
    \item{\code{D_GENE_and_allele}}{a character vector}
    \item{\code{D_REGION_reading_frame}}{a numeric vector}
    \item{\code{CDR1_IMGT_length}}{a numeric vector}
    \item{\code{CDR2_IMGT_length}}{a character vector}
    \item{\code{CDR3_IMGT_length}}{a character vector}
    \item{\code{CDR_IMGT_lengths}}{a character vector}
    \item{\code{FR_IMGT_lengths}}{a character vector}
    \item{\code{AA_JUNCTION}}{a character vector}
    \item{\code{JUNCTION_frame}}{a character vector}
    \item{\code{Orientation}}{a character vector}
    \item{\code{Functionality_comment}}{a character vector}
    \item{\code{V_REGION_potential_ins/del}}{a character vector}
    \item{\code{J_GENE_and_allele_comment}}{a character vector}
    \item{\code{V_REGION_insertions}}{a character vector}
    \item{\code{V_REGION_deletions}}{a character vector}
    \item{\code{Sequence}}{a character vector}
  }
}

\references{
Alamyar, E. et al., IMGT/HighV-QUEST: A High-Throughput System and Web Portal for the Analysis of Rearranged Nucleotide Sequences of Antigen Receptors, JOBIM 2010, Paper 63 (2010).

Brochet, X. et al., Nucl. Acids Res. 36, W503-508 (2008). PMID: 18503082

IMGT/HighV-QUEST 1_Summary file description: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#Esummary}
}
\examples{
data(summarytab)
str(summarytab)
}
\keyword{datasets}
