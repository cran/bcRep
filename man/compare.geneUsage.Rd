\name{compare.geneUsage}
\alias{compare.geneUsage}
\alias{plotCompareGeneUsage}
\title{
Compare gene usage of different samples
}
\description{
This function compares the gene usage of different samples (see details!). Analysis can be done for subgroups, genes and alleles (see Details). Values can be returned as relative or absolute abundance. Parallel processing is possible.
}
\usage{
compare.geneUsage(gene.list = NULL, level = c("subgroup", "gene", "allele"), 
     abundance = c("relative", "absolute"), names = NULL, nrCores = 1)

plotCompareGeneUsage(comp.tab = NULL, color = c("gray97", "darkblue"), 
     title = NULL, PDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.list}{
     A list containing vectors of genes of each sample
}
  \item{level}{
	Gene level used for gene usage analysis: subgroup, gene, allele
}
  \item{abundance}{
	Shall relative or absolute values be returned? (default: relative)
}
  \item{names}{
     A vector containing names for the samples (default: like Sample1, Sample2, ...)
}
  \item{nrCores}{
     Number of cores used for parallel processing
}
  \item{comp.tab}{
     Output tab from \code{compare.geneUsage()}
}
  \item{color}{
     Colors used for heatmap (default: gray and darkblue)
}
  \item{title}{
     Title of plot
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
Gene usage analysis will be done for each sample. Vector of genes will be analyzed for one of the levels \code{subgroup} (e.g. IGHV1), \code{gene} (e.g. IGHV1-1) or \code{allele} (e.g. IGHV1-1*2). Proportions (\code{abundance = "relative"}) are always based on the number of all alleles found in list: the number of the subgroup/gene/allele is divided by the number of all alleles mentioned for all sequences (in the case there are more alleles/genes mentioned for one sequence).


The \code{PDF} character string should be only the project name (without ".pdf"). 

A figure called "PDF"_Comparison_Gene-usage.pdf will be saved to the working directory.
}
\value{
Output is a data frame containing absolue or relative values of gene usage of each sample.
}
\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
\code{\link{geneUsage}}, \code{\link{compare.geneUsage}}, \code{\link{plotCompareGeneUsage}}
}
\examples{
data(aaseqtab)
data(aaseqtab2)

Vgenes.comp<-compare.geneUsage(gene.list = list(aaseqtab$V_GENE_and_allele, 
     aaseqtab2$V_GENE_and_allele), level = "subgroup", abundance = "relative", 
     names = c("IndA", "IndB"), nrCores = 1)
\dontrun{
plotCompareGeneUsage(comp.tab = Vgenes.comp, color = c("gray97", "darkblue"), PDF = "Example")
}
}
