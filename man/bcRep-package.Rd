\name{bcRep-package}
\alias{bcRep-package}
\alias{bcRep}
\docType{package}
\title{
Advanced Analysis of B Cell Receptor Repertoire Data

}
\description{
This package helps to analyze IMGT/HighV-QUEST output data, in more detail. It functions well with B cells, but can also be used for T cell data, in some cases.
Using this package IMGT/HighV-QUEST output files can be readed and sequences and clones studied. In special their functionality, junction frames, gene usage and mutations. Functions to analyze clones, but also to compare sequences and clones, are provided. Further distances based on sequence or gene usage data can be calculated and multidimensional scaling can be performed.
}
\details{
\tabular{ll}{
Package: \tab bcRep\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-09\cr
License: GPL-2
}
For many of the functions output data from IMGT/HighV-QUEST is required. In this case the particularly input file is given in the corresponding help file and \code{readIMGT()} can be used to read these files. IMGT files from different folders can be combined using \code{combineIMGT()}. For sequence analysis functions like \code{sequences.functionality()} or \code{sequences.junctionFrame()} exist, which give an overview about functionality or junction frame usage of sequences. Further gene usage (\code{geneUsage()}) and gene/gene combinations (\code{sequences.geneComb()}) can be analyzed. The function \code{sequences.mutation()} returns an overview of all mutations, replacement and silent mutations, as well as the R/S ratio in several regions like V, FR1-3, CDR1-2. IMGT/HighV-QUEST output can also be used to analyze clones (\code{clones()}). Therefore criteria can be changed (CDR3 identity, V and J gene usage) and results of different samples be compared (\code{clones.shared()}). Further more amino acid distributions, as well as Gini index, richness and diversity of sequences, can be studied (\code{AADistribution()}, \code{trueDiversity()}, \code{clones.GiniIndex()}). Principal coordinate analysis on sequences, as well as on gene usage data can be perfomed, using different distances (\code{sequences.distance()}, \code{geneUsage.distance()}, \code{dist.PCoA}). Special plot function for most of the methods are provided.
}
\author{
Julia Bischof

Maintainer: Julia Bischof <Julia.Bischof@uksh.de>
}
\references{
Alamyar, E. et al., IMGT/HighV-QUEST: A High-Throughput System and Web Portal for the Analysis of Rearranged Nucleotide Sequences of Antigen Receptors, JOBIM 2010, Paper 63 (2010). Website: \url{http://www.imgt.org/}

Brochet, X. et al., Nucl. Acids Res. 36, W503-508 (2008). PMID: 18503082

IMGT/V-QUEST Documentation: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#output3}

IMGT Repertoire (IG and TR): \url{http://www.imgt.org/IMGTrepertoire/LocusGenes/}

Lou Jost: Entropy and diversity; OIKOS 113:2 (2006)
}

\keyword{ package }

\examples{
\dontrun{
data(summarytab)
data(aaseqtab)
data(aaseqtab2)
data(mutationtab)
data(clones.ind)
data(clones.allind)
data(vgenes) 

## Combine IMGT/HighV-QUEST folders and read data
combineIMGT(folders = c("pathTo/IMGT1a", "pathTo/IMGT1b", "pathTo/IMGT1c"), 
name = "NewProject)
tab<-readIMGT("PathTo/file.txt",filterNoResults=TRUE)


## Get information about functionality and filter for functional sequences
functionality<-sequences.functionality(data = summarytab$Functionality)
ProductiveSequences<-sequences.getProductives(summarytab)

## Gene usage
Vsubgroup.usage<-geneUsage(genes = clones.ind$V_gene, 
     functionality = clones.ind$Functionality_all_sequences, level = "subgroup", 
     abundance="relative")

Vgenes.comp<-compare.geneUsage(gene.list = list(aaseqtab$V_GENE_and_allele, 
     aaseqtab2$V_GENE_and_allele), level = "subgroup", abundance = "relative", 
     names = c("IndA", "IndB"), nrCores = 1)
plotCompareGeneUsage(comp.tab = Vgenes.comp, color = c("gray97", "darkblue"), 
     PDF = "Example")


## Gene/gene combinations
VDcomb.tab<-sequences.geneComb(family1 = summarytab$V_GENE_and_allele, 
     family2 = summarytab$D_GENE_and_allele, level = "subgroup", 
     abundance = "relative")
plotGeneComb(geneComb.tab=VDcomb.tab, color="red", withNA=FALSE,PDF="test")


## Mutation analysis
mutation.V<-sequences.mutation(mutationtab = mutationtab, summarytab = summarytab, 
     sequence = "V")
mutation.CDR1<-sequences.mutation(mutationtab = mutationtab, sequence = "CDR1", 
     functionality = TRUE, junctionFr = TRUE)

## Defining, Filtering and Plotting Clone features
clones.tab<-clones(aaseqtab=aaseqtab,summarytab=summarytab, identity=0.85, useJ=TRUE, 
     dispCDR3aa=TRUE, dispFunctionality.ratio=TRUE, dispFunctionality.list=TRUE)
plotClonesCDR3Length(CDR3Length = clones.ind$CDR3_length_AA, 
     functionality = clones.ind$Functionality_all_sequences,  
     color="gray",abundance="relative", PDF="test")
clones.func<-clones.filterFunctionality(clones.tab = clones.ind, 
     filter = "productive")   

## Find shared clones between individuals
sharedclones<-clones.shared(clones.tab = clones.allind, identity = 0.85, useJ = TRUE, 
     dispD = TRUE, dispCDR3aa = TRUE)
sharedclones.summary<-clones.shared.summary(shared.tab = sharedclones)

## True diversity
trueDiv<-trueDiversity(sequences = aaseqtab$CDR3_IMGT, order = 1)
plotTrueDiversity(trueDiversity.tab=trueDiv,color="red",PDF="test")

trueDiv.comp<-compare.trueDversity(sequence.list = list(aaseqtab$CDR3_IMGT, 
     aaseqtab2$CDR3_IMGT), names = c("IndA", "IndB"), order = 1, nrCores = 1)
plotCompareTrueDiversity(comp.tab = trueDiv.comp, PDF = "Example")


## Gini index
gini<-gini<-clones.giniIndex(clone.size=clones.ind$total_number_of_sequences)

## Dissmilarity/distance indices of gene usage and sequence data

distGeneUsage<-geneUsage.distance(geneUsage.tab = Vgenes, method = "bc")
distSequence<-sequences.distance(sequences = clones.ind$unique_CDR3_sequences_AA, 
     method = "levenshtein", divLength=TRUE)

## Principal coordinate analysis of distance matrices + visualization
distpcoa<-dist.PCoA(dist.tab = distGeneUsage, correction = "none")
     # 'groups' data.frame for plot function: in the case, there are no groups:
groups.vec<-unlist(apply(data.frame(clones.ind$unique_CDR3_sequences_AA),1,
          function(x){strsplit(x,split=", ")[[1]]}))
groups.vec<-cbind(groups.vec, 1)
plotDistPCoA(pcoa.tab = distpcoa, groups=groups.vec, axes = c(1,2), 
     plotCorrection = FALSE, title = NULL, plotLegend=TRUE, PDF = "TEST")    

}
}
