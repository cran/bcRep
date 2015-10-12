\name{bcRep-package}
\alias{bcRep-package}
\alias{bcRep}
\docType{package}
\title{
Advanced Analysis of B Cell Receptor Repertoire Data

}
\description{
This package helps to analyze IMGT/HighV-QUEST output, in more detail. It functions well with B cells, but can also be used for T cell data, in some cases.
Using this package you can read IMGT/HighV-QUEST output files and study sequences and clones. In special their functionality, junction frames, gene usage and mutations. Functions to analyze clones out of the IMGT/HighV-QUEST output, but also to compare shared clones, are provided.
}
\details{
\tabular{ll}{
Package: \tab bcRep\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-09\cr
License: GPL-2
}
For many of the functions you need output data from IMGT/HighV-QUEST. In this case the required input file is given in help and you can use \code{readIMGT()} to read these files. For sequence analysis you can use functions like \code{sequences.functionality()} or \code{sequences.junctionFrame()}, which give an overview about functionality or junction frame usage of sequences. Further gene usage (\code{geneUsage()}) and gene/gene combinations (\code{sequences.geneComb()}) can be analyzed. The function \code{sequences.mutation()} returns an overview of all mutations, replacement and silent mutations, as well as the R/S ratio in several regions like V, FR1-3, CDR1-2. IMGT/HighV-QUEST output can also be used to analyze clones (\code{clones()}). Therefore criteria can be changed and results of for example different individuals be compared (\code{clones.shared()}). Further more amino acid distributions, as well as richness and diversity of sequences, like CDR3, can be studied (\code{AADistribution()}, \code{trueDiversity()}). Special plot function for most of the methods are provided.
}
\author{
Julia Bischof

Maintainer: Julia Bischof <Julia.Bischof@uksh.de>
}
\references{
Alamyar, E. et al., IMGT/HighV-QUEST: A High-Throughput System and Web Portal for the Analysis of Rearranged Nucleotide Sequences of Antigen Receptors, JOBIM 2010, Paper 63 (2010).

Brochet, X. et al., Nucl. Acids Res. 36, W503-508 (2008). PMID: 18503082

IMGT/V-QUEST Documentation: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#output3}

IMGT Repertoire (IG and TR): \url{http://www.imgt.org/IMGTrepertoire/LocusGenes/}

Lou Jost: Entropy and diversity; OIKOS 113:2 (2006)
}

\keyword{ package }

\examples{
data(summarytab)
data(aaseqtab)
data(mutationtab)
data(clones.ind)
data(clones.allind)

## Read IMGT/HighV-QUEST data
\dontrun{tab<-readIMGT("PathTo/file.txt",filterNoResults=TRUE)}

## Get information about functionality and filter for functional sequences
functionality<-sequences.functionality(data = summarytab$Functionality)
ProductiveSequences<-sequences.getProductives(summarytab)

## Gene usage
Vsubgroup.usage<-geneUsage(geneUsage = clones.ind$V_gene, 
     functionality = clones.ind$Functionality_all_sequences, level = "subgroup", 
     abundance="relative")

## Gene/gene combinations
\dontrun{VDcomb.tab<-sequences.geneComb(family1 = summarytab$V_GENE_and_allele, 
     family2 = summarytab$D_GENE_and_allele, level = "subgroup", abundance = "relative")
plotGeneComb(geneComb.tab=VDcomb.tab, color="red", withNA=FALSE,PDF="test")}

## Mutation analysis
mutation.V<-sequences.mutation(mutationtab = mutationtab, summarytab = summarytab, 
     sequence = "V")
mutation.CDR1<-sequences.mutation(mutationtab = mutationtab, sequence = "CDR1", 
     functionality = TRUE, junctionFr = TRUE)

## Define clones and CDR3 length analysis of clones
\dontrun{clones.tab<-clones(aaseqtab=aaseqtab,summarytab=summarytab, identity=0.85, useJ=TRUE, 
     dispCDR3aa=TRUE, dispFunctionality.ratio=TRUE, dispFunctionality.list=TRUE)
plotClonesCDR3Length(CDR3Length = clones.ind$CDR3_length_AA, 
     functionality = clones.ind$Functionality_all_sequences,  
     color="gray",abundance="relative", PDF="test")}

## Find shared clones between individuals
sharedclones<-clones.shared(clones.tab = clones.allind, identity = 0.85, useJ = TRUE, 
     dispD = TRUE, dispCDR3aa = TRUE)
sharedclones.summary<-clones.shared.summary(shared.tab = sharedclones)

## True diversity
trueDiv<-trueDiversity(sequences = aaseqtab$CDR3_IMGT, order = 1)
\dontrun{plotTrueDiversity(trueDiversity.tab=trueDiv,color="red",PDF="test")}

}
