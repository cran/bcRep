\name{clones}
\alias{clones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grouping sequences into clones
}
\description{
This function uses IMGT/HighV-QUEST output files to define B cell clones. Therefore criteria using amino acid CDR3 sequences, V genes and J genes (optional) are used. A treshold for CDR3 identity can be given. Parallel processing is possible.
}
\usage{
clones(aaseqtab = NULL, summarytab = NULL, ntseqtab = NULL, identity = 0.85, 
     useJ = TRUE,dispD = FALSE, dispSeqID = FALSE, dispCDR3aa = FALSE, 
     dispCDR3nt = FALSE, dispJunctionFr.ratio = FALSE, 
     dispJunctionFr.list = FALSE, dispFunctionality.ratio = FALSE, 
     dispFunctionality.list = FALSE, dispTotalSeq = FALSE, nrCores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aaseqtab}{
	IMGT/HighV-QUEST output, file 5_AA-sequences(...).txt
}
  \item{summarytab}{
	IMGT/HighV-QUEST output, file 1_Summary(...).txt
}
  \item{ntseqtab}{
	IMGT/HighV-QUEST output, file 3_Nt-sequences(...).txt (optional)
}
  \item{identity}{
	Treshold of CDR3 identity. A value between 0 and 1.
}
  \item{useJ}{
	Shall J gene be included into analysis? default: TRUE
}
  \item{dispD}{
	Shall D genes and alleles be returned? default: FALSE
}
  \item{dispSeqID}{
	Shall sequence ID's be returned? default: FALSE
}
  \item{dispCDR3aa}{
	Shall amino acid CDR3 sequences be returned? default: FALSE
}
  \item{dispCDR3nt}{
	Shall nucleotide amino acid sequences be returned? default: FALSE
}
  \item{dispJunctionFr.ratio}{
	Shall ratios of in-frame, out-of-frame and unknown junctions be returned? default: FALSE
}
  \item{dispJunctionFr.list}{
	Shall a list of all junction frames be returned? default: FALSE
}
  \item{dispFunctionality.ratio}{
	Shall ratios of productive, unproductive and unknown functionality sequences be returned? default: FALSE
}
  \item{dispFunctionality.list}{
	Shall a list of all functionalities be returned? default: FALSE
}
  \item{dispTotalSeq}{
	Shall all total nucleotide sequences be returned? default: FALSE
}
  \item{nrCores}{
	Number of cores used for parallel processing (default: 1)
}
}
\details{
This function uses IMGT/HighV-QUEST output to define clones. Therefore amino acid CDR3 sequences, V genes and J genes (optional) are used. 
Criteria for clone groups are 1) same CDR3 length, 2) CDR3 identity of a given treshold, 3) same V gene and 4) same J gene (optional). 
A treshold for CDR3 identity has to be between 0 and 1. A cutoff of 0.85 means CDR3 identity of 85\%. For example for a CDR3 length of 15 amino acids 85\% identity would mean that at least 11 of 15 positions have to be identical (0.85*15 = 10.75; values are rounded up).

\code{useJ=T} includes also the criteria of same J genes for clone defintion.

}
\value{
Output of clones() is a data frame containing
  \item{unique_CDR3_sequences_[AA]}{unique CDR3 sequences belonging to this clone}
  \item{CDR3_length_AA}{CDR3 length in amino acids}
  \item{number_of_unique_sequences}{number of unique CDR3 sequences belonging to this clone}
  \item{total_number_of_sequences}{number of all sequences belonging to this clone (one sequence can appear several times)}
  \item{sequence_count_per_CDR3}{sequence count for each of the unique CDR3 sequences}
  \item{V_gene}{V gene belonging to this clone}
  \item{V_gene_and_allele}{original IMGT V gene nomenclature}
  \item{J_gene}{J gene(s) belonging to this clone  (if useJ=F, there can be several J genes)}
  \item{J_gene_and_allele}{original IMGT J gene nomenclature}
  \item{optional arguments}{D_gene;_all_CDR3_sequences_AA; all_CDR3_sequences_nt; 
  Funct_all_sequences; 
  Funct_productive/unproductive/unknown sequences; 
  Junction_frame_all_sequences; 
  JF_in-frame/out-of-frame/unknown sequences; 
  Sequence_IDs; Total_sequences_nt}
}
\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
 \code{\link{clones.CDR3Length}}, \code{\link{plotClonesCDR3Length}}, \code{\link{plotClonesCopyNumber}}, \code{\link{geneUsage}}, 
 
 \code{\link{plotGeneUsage}}, \code{\link{clones.shared}}
}
\examples{
data(summarytab)
data(aaseqtab)

clones.tab<-clones(aaseqtab=aaseqtab,summarytab=summarytab, identity=0.85, useJ=TRUE, 
     dispCDR3aa=TRUE, dispFunctionality.ratio=TRUE, dispFunctionality.list=TRUE)
}

