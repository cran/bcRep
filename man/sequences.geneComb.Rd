\name{sequences.geneComb}
\alias{sequences.geneComb}
\alias{plotGeneComb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene/gene combinations
}
\description{
This function calculates the abundance of two paired gene families. Analysis can be done for subgroup, gene and allele level (see Details).

\code{sequences.geneComb} returns a data frame containing all possible combinations of gene families 1 (rows) and 2 (columns) in relative or absolute values.

\code{plotGeneComb} visualizes relative abundances of gene/gene combinations in a heatmap. Absolute values will be converted to relative abundances.
}
\usage{
sequences.geneComb(family1 = NULL, family2 = NULL, 
     level = c("subgroup","gene","allele"), abundance = c("relative", "absolute"),
     nrCores=1)
     
plotGeneComb(geneComb.tab=NULL,color=c("gray97","darkblue"),withNA=TRUE, 
     title=NULL, PDF=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family1}{
     Vector including gene family 1, should be matched to order of family2 (see Details)
}
  \item{family2}{
     Vector including gene family 2, should be matched to order of family1 (see Details)
}
  \item{level}{
     Gene family level: subgroup, gene or allele (see Details)
}
  \item{abundance}{
     How values shall be returned: relative or absolute abundance.
}
  \item{nrCores}{
     Number of cores used for parallel processing (default: 1)
}
  \item{geneComb.tab}{
     Output data frame from \code{sequences.genecomb()}
}
  \item{color}{
     colors of heatmap (default: gray and blue)
}
  \item{withNA}{
     Shall combination without information for one of both families be included? (default: TRUE)
}
  \item{title}{
     title of plot (NULL: "family1 & family2 gene combinations (zoom: xx)"; turn of with "")
}
\item{PDF}{
     PDF project name (see Details)
}
\item{...}{

}
}
\details{
Input are vectors of IMGT/HighV-QUEST output file 1_Summary(...).txt. Columns like "V-GENE and allele", "D-GENE and allele" or "J-GENE and allele" can be used as input, but need to have the same order! Levels can be \code{subgroup} (e.g. IGHV1), \code{genes} (e.g. IGHV1-1) or \code{alleles} (e.g. IGHV1-1*2).

The \code{PDF} character string should be only the project name (without ".pdf"). 

A figure called "PDF"_"family1"-"family2"-combinations.pdf will be saved to the working directory.
}
\value{
Output is a data frame containing gene/gene combinations, with family 1 as rows and family 2 as columns. Either in absolute or relative abundances.
}
\references{
IMGT Repertoire (IG and TR): \url{http://www.imgt.org/IMGTrepertoire/LocusGenes/}
}
\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
     \code{\link{sequences.geneComb}}, \code{\link{plotGeneComb}}, \code{\link{geneUsage}}
}
\examples{
data(summarytab)
VDcomb.tab<-sequences.geneComb(family1 = summarytab$V_GENE_and_allele, 
     family2 = summarytab$D_GENE_and_allele, level = "subgroup", abundance = "relative")
     
plotGeneComb(geneComb.tab=VDcomb.tab, withNA=FALSE)
}

