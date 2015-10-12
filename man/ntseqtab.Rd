\name{ntseqtab}
\alias{ntseqtab}
\docType{data}
\title{
Summary of nucleotide sequences of B cells of one individual
}
\description{
Data frame represents an output from IMGT/HighV-QUEST, file 3_Nt-sequences(...).txt
}
\usage{data("ntseqtab")}
\format{
  A data frame with 3000 observations on the following 114 variables.
  \describe{
    \item{\code{Sequence_number}}{a numeric vector}
    \item{\code{Sequence_ID}}{a character vector}
    \item{\code{Functionality}}{a character vector}
    \item{\code{V_GENE_and_allele}}{a character vector}
    \item{\code{J_GENE_and_allele}}{a character vector}
    \item{\code{D_GENE_and_allele}}{a character vector}
    \item{\code{V_D_J_REGION}}{a character vector}
    \item{\code{V_J_REGION}}{a character vector}
    \item{\code{V_REGION}}{a character vector}
    \item{\code{FR1_IMGT}}{a character vector}
    \item{\code{CDR1_IMGT}}{a character vector}
    \item{\code{FR2_IMGT}}{a character vector}
    \item{\code{CDR2_IMGT}}{a character vector}
    \item{\code{FR3_IMGT}}{a character vector}
    \item{\code{CDR3_IMGT}}{a character vector}
    \item{\code{JUNCTION}}{a character vector}
    \item{\code{3V_REGION}}{a character vector}
    \item{\code{(N_D)_J_REGION}}{a character vector}
    \item{\code{(N_D)_REGION}}{a character vector}
    \item{\code{P3V}}{a character vector}
    \item{\code{N_REGION}}{a character vector}
    \item{\code{N1_REGION}}{a character vector}
    \item{\code{P5D}}{a character vector}
    \item{\code{D_REGION}}{a character vector}
    \item{\code{P3D}}{a character vector}
    \item{\code{P5D1}}{a logical vector}
    \item{\code{D1_REGION}}{a character vector}
    \item{\code{P3D1}}{a character vector}
    \item{\code{N2_REGION}}{a character vector}
    \item{\code{P5D2}}{a character vector}
    \item{\code{D2_REGION}}{a character vector}
    \item{\code{P3D2}}{a character vector}
    \item{\code{N3_REGION}}{a character vector}
    \item{\code{P5D3}}{a character vector}
    \item{\code{D3_REGION}}{a character vector}
    \item{\code{P3D3}}{a character vector}
    \item{\code{N4_REGION}}{a character vector}
    \item{\code{P5J}}{a character vector}
    \item{\code{5J_REGION}}{a character vector}
    \item{\code{D_J_REGION}}{a character vector}
    \item{\code{J_REGION}}{a character vector}
    \item{\code{FR4_IMGT}}{a character vector}
    \item{\code{V_D_J_REGION_start}}{a numeric vector}
    \item{\code{V_D_J_REGION_end}}{a numeric vector}
    \item{\code{V_J_REGION_start}}{a numeric vector}
    \item{\code{V_J_REGION_end}}{a logical vector}
    \item{\code{V_REGION_start}}{a numeric vector}
    \item{\code{V_REGION_end}}{a numeric vector}
    \item{\code{FR1_IMGT_start}}{a numeric vector}
    \item{\code{FR1_IMGT_end}}{a numeric vector}
    \item{\code{CDR1_IMGT_start}}{a numeric vector}
    \item{\code{CDR1_IMGT_end}}{a numeric vector}
    \item{\code{FR2_IMGT_start}}{a numeric vector}
    \item{\code{FR2_IMGT_end}}{a numeric vector}
    \item{\code{CDR2_IMGT_start}}{a numeric vector}
    \item{\code{CDR2_IMGT_end}}{a numeric vector}
    \item{\code{FR3_IMGT_start}}{a numeric vector}
    \item{\code{FR3_IMGT_end}}{a numeric vector}
    \item{\code{CDR3_IMGT_start}}{a numeric vector}
    \item{\code{CDR3_IMGT_end}}{a numeric vector}
    \item{\code{JUNCTION_start}}{a numeric vector}
    \item{\code{JUNCTION_end}}{a numeric vector}
    \item{\code{3V_REGION_start}}{a numeric vector}
    \item{\code{3V_REGION_end}}{a numeric vector}
    \item{\code{(N_D)_J_REGION_start}}{a numeric vector}
    \item{\code{(N_D)_J_REGION_end}}{a numeric vector}
    \item{\code{(N_D)_REGION_start}}{a numeric vector}
    \item{\code{(N_D)_REGION_end}}{a numeric vector}
    \item{\code{P3V_start}}{a numeric vector}
    \item{\code{P3V_end}}{a numeric vector}
    \item{\code{N_REGION_start}}{a numeric vector}
    \item{\code{N_REGION_end}}{a numeric vector}
    \item{\code{N1_REGION_start}}{a numeric vector}
    \item{\code{N1_REGION_end}}{a numeric vector}
    \item{\code{P5D_start}}{a numeric vector}
    \item{\code{P5D_end}}{a numeric vector}
    \item{\code{D_REGION_start}}{a numeric vector}
    \item{\code{D_REGION_end}}{a numeric vector}
    \item{\code{P3D_start}}{a numeric vector}
    \item{\code{P3D_end}}{a numeric vector}
    \item{\code{P5D1_start}}{a numeric vector}
    \item{\code{P5D1_end}}{a numeric vector}
    \item{\code{D1_REGION_start}}{a numeric vector}
    \item{\code{D1_REGION_end}}{a numeric vector}
    \item{\code{P3D1_start}}{a numeric vector}
    \item{\code{P3D1_end}}{a numeric vector}
    \item{\code{N2_REGION_start}}{a numeric vector}
    \item{\code{N2_REGION_end}}{a numeric vector}
    \item{\code{P5D2_start}}{anumeric vector}
    \item{\code{P5D2_end}}{a numeric vector}
    \item{\code{D2_REGION_start}}{a numeric vector}
    \item{\code{D2_REGION_end}}{a numeric vector}
    \item{\code{P3D2_start}}{a numeric vector}
    \item{\code{P3D2_end}}{a numeric vector}
    \item{\code{N3_REGION_start}}{a numeric vector}
    \item{\code{N3_REGION_end}}{a numeric vector}
    \item{\code{P5D3_start}}{a numeric vector}
    \item{\code{P5D3_end}}{a numeric vector}
    \item{\code{D3_REGION_start}}{a numeric vector}
    \item{\code{D3_REGION_end}}{a numeric vector}
    \item{\code{P3D3_start}}{a numeric vector}
    \item{\code{P3D3_end}}{a numeric vector}
    \item{\code{N4_REGION_start}}{a numeric vector}
    \item{\code{N4_REGION_end}}{a numeric vector}
    \item{\code{P5J_start}}{a numeric vector}
    \item{\code{P5J_end}}{a numeric vector}
    \item{\code{5J_REGION_start}}{a numeric vector}
    \item{\code{5J_REGION_end}}{a numeric vector}
    \item{\code{D_J_REGION_start}}{a numeric vector}
    \item{\code{D_J_REGION_end}}{a numeric vector}
    \item{\code{J_REGION_start}}{a numeric vector}
    \item{\code{J_REGION_end}}{a numeric vector}
    \item{\code{FR4_IMGT_start}}{a numeric vector}
    \item{\code{FR4_IMGT_end}}{a numeric vector}
  }
}

\references{
Alamyar, E. et al., IMGT/HighV-QUEST: A High-Throughput System and Web Portal for the Analysis of Rearranged Nucleotide Sequences of Antigen Receptors, JOBIM 2010, Paper 63 (2010).

Brochet, X. et al., Nucl. Acids Res. 36, W503-508 (2008). PMID: 18503082

IMGT/HighV-QUEST 3_Nt-sequences file description: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#Ent}
}
\examples{
data(ntseqtab)
str(ntseqtab)
}
\keyword{datasets}
