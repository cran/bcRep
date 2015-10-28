\name{clones.shared}
\alias{clones.shared}
\alias{clones.shared.summary}
\title{
Shared clones between individuals
}
\description{
\code{clones.shared} compares clones of different individuals to find shared clones with identical or similar CDR3 sequences. Criteria for same/similar clones are 1) same CDR3 length, 2) CDR3 identity of a given treshold, 3) same V gene, 4) same J gene (optional). Parallel processing is possible.

\code{clones.shared.summary} provides information about the number of shared clones between two or more individuals.
}
\usage{
clones.shared(clones.tab = NULL, identity = 0.85, useJ = TRUE, dispD = TRUE, 
     dispCDR3aa = TRUE, dispCDR3nt = FALSE, dispFunctionality.list = FALSE, 
     dispFunctionality.ratio = FALSE, dispJunctionFr.list = FALSE, 
     dispJunctionFr.ratio = FALSE, dispTotalSeq = FALSE, nrCores=1)

clones.shared.summary(shared.tab = NULL, clones.tab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clones.tab}{
     A dataframe, which includes all individual clones. Therefore data frames, as output from \code{clones()} has to be combined, via \code{rbind()}. The first column needs to include individual identifiers. Example: \code{data(ClonesAllInd)}
}
  \item{identity}{
     A Value between 0 and 1, indicating proportion of identity of sequences.
}
  \item{useJ}{
     TRUE: J gene shall be included as criteria (default: TRUE)
}
  \item{dispD}{
     TRUE: return D gene and allele in output data frame (only, if they are included in input; default: TRUE)
}
  \item{dispCDR3aa}{
     TRUE: return CDR3 amino acid sequences in output data frame (only, if they are included in input; default: TRUE)
}
  \item{dispCDR3nt}{
     TRUE: return CDR3 nucleotide sequences in output data frame (only, if they are included in input; default: FALSE)
}
  \item{dispFunctionality.list}{
     TRUE: return functionality list in output data frame (only, if they are included in input; default: FALSE)
}
  \item{dispFunctionality.ratio}{
     TRUE: return functionality ratios for prosuctive and unproductive sequences in output data frame (only, if they are included in input; default: FALSE)
}
  \item{dispJunctionFr.list}{
     TRUE: return junction frame list in output data frame (only, if they are included in input; default: FALSE)
}
  \item{dispJunctionFr.ratio}{
     TRUE: return junction frame ratios for in-frame and out-of-frame sequences output data frame (only, if they are included in input; default: FALSE)
}
  \item{dispTotalSeq}{
     TRUE: return total nucleotide sequences output data frame (only, if they are included in input; default: FALSE)
}
  \item{nrCores}{
	Number of cores used for parallel processing (default: 1)
}
  \item{shared.tab}{
     Output from function \code{clones.shared()}. Either whole data frame or individual column.
}
}
\details{
\code{clones.shared}: Criteria for clone groups are 1) same CDR3 length, 2) CDR3 identity of a given treshold, 3) same V gene and 4) same J gene (optional). Same or similar CDR3 sequence has to be shared between at least two individuals.
A treshold for CDR3 identity has to be between 0 and 1. A cutoff of 0.85 means CDR3 identity of 85\%. For example for a CDR3 length of 15 amino acids 85\% identity would mean that at least 11 of 15 positions have to be identical (0.85*15 = 10.75; values are rounded up).

\code{useJ=T} includes also the criteria of same J genes for clone defintion.

\code{clones.shared.summary} summarizes information about shared clones. If \code{clonestab.individual} is also provided, number of clones, appearing in only one of these individuals is also returned.
}
\value{
Output of \code{clones.shared} is a data frame including
  \item{number_individuals}{The number of individuals, that share this clones}  
  \item{individuals}{Individual ID's, seperated by ";"}
  \item{CDR3_length_AA}{length of CDR3 amino acid sequence}
  \item{shared_CDR3}{shared (100\% identity) or similar (<100\% identity) CDR3 sequence}
  \item{number_shared_CDR3}{Number of shared CDR3 sequences}  
  \item{unique_CDR3_sequences_AA_per_individual}{CDR3 sequences per clone, individuals are seperated by ";"}  
  \item{sequence_count_per_CDR3}{Quantity how often clones of "unique all CDR3 sequences [AA] per individual" appaer, individuals are seperated by ";"}  
  \item{V_gene}{V gene}  
  \item{J_gene}{J gene(s)}  
  \item{optional output}{D_gene; all_CDR3_sequences_AA; al_CDR3_sequences_nt; Funct_all_sequences; Funct_productive/unproductive; Junction_frame_all_sequences; JF_in-frame/out-of-frame; Total_sequences_nt (individuals are seperated by ";")}

Output of \code{clones.shared.summary} is a data frame containing all possible groups and the quantity of appearance.
}

\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
\code{\link{clones}}
}
\examples{
data(clones.allind)
sharedclones<-clones.shared(clones.tab = clones.allind, identity = 0.85, useJ = TRUE, 
     dispD = TRUE, dispCDR3aa = TRUE)
sharedclones.summary<-clones.shared.summary(shared.tab = sharedclones)
}

