\name{aaseqtab}
\alias{aaseqtab}
\docType{data}
\title{
Summary of amino acid sequences of B cells of one individual
}
\description{
Data frame represents an output from IMGT/HighV-QUEST, file 5_AA-sequences(...).txt
}
\usage{data("aaseqtab")}
\format{
  A data frame with 3000 observations on the following 18 variables.
  \describe{
    \item{\code{Sequence_number}}{a numeric vector}
    \item{\code{Sequence_ID}}{a character vector}
    \item{\code{Functionality}}{a character vector}
    \item{\code{V_GENE_and_allele}}{a character vector}
    \item{\code{J_GENE_and_allele}}{a character vector}
    \item{\code{D_GENE_and_allele}}{a character vector}
    \item{\code{V_D_J_REGION}}{a character vector}
    \item{\code{V_J_REGION}}{a logical vector}
    \item{\code{V_REGION}}{a character vector}
    \item{\code{FR1_IMGT}}{a character vector}
    \item{\code{CDR1_IMGT}}{a character vector}
    \item{\code{FR2_IMGT}}{a character vector}
    \item{\code{CDR2_IMGT}}{a character vector}
    \item{\code{FR3_IMGT}}{a character vector}
    \item{\code{CDR3_IMGT}}{a character vector}
    \item{\code{JUNCTION}}{a character vector}
    \item{\code{J_REGION}}{a character vector}
    \item{\code{FR4_IMGT}}{a character vector}
  }
}
\references{
Alamyar, E. et al., IMGT/HighV-QUEST: A High-Throughput System and Web Portal for the Analysis of Rearranged Nucleotide Sequences of Antigen Receptors, JOBIM 2010, Paper 63 (2010).

Brochet, X. et al., Nucl. Acids Res. 36, W503-508 (2008). PMID: 18503082

IMGT/HighV-QUEST 5_AA-sequences file description: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#EAA}
}
\examples{
data(aaseqtab)
str(aaseqtab)
}
\keyword{datasets}
