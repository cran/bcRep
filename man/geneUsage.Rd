\name{geneUsage}
\alias{geneUsage}
\alias{plotGeneUsage}
\title{
Gene usage statistics
}
\description{
This function gives information about the gene usage distribution (see detail!). Results can returned as relative or absolute values and be visualized as a barplot, using \code{plotGeneUsage}. Further gene usage can be analyzed in relation to functionality or junction frame usage.
}
\usage{
geneUsage(genes = NULL, level = c("subgroup", "gene", "allele"), 
     functionality = NULL, junctionFr = NULL, 
     abundance=c("relative","absolute"), ...)

plotGeneUsage(geneUsage.tab=NULL,plotFunctionality=FALSE,plotJunctionFr=FALSE,
     color=c("orange","darkblue","gray"), title=NULL, PDF=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{genes}{
	Vector containing gene information (IMGT nomenclature; see Details)
}
\item{level}{
	Gene level used for gene usage analysis: subgroup, gene, allele
}
\item{functionality}{
	Vector containing functionality information
}
\item{junctionFr}{
	Vector containing functionality information
}
\item{abundance}{
	Shall relative or absolute values be returned? (default: relative)
}
\item{geneUsage.tab}{
	Output of \code{geneUsage()}
}
\item{plotFunctionality}{
	Shall gene usage vs. functionality be plotted? (default: FALSE)
}
\item{plotJunctionFr}{
     Shall gene usage vs. junction frame be plotted? (default: FALSE)
}
\item{color}{
	color used for plots (default: c("orange","darkblue","gray"))
}

\item{title}{
	title of the plot (optional)
}
\item{PDF}{
     PDF project name (see Details)
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Vector of genes will be analyzed for one of the levels \code{subgroup} (e.g. IGHV1), \code{gene} (e.g. IGHV1-1) or \code{allele} (e.g. IGHV1-1*2). Proportions (\code{abundance = "relative"}) are always based on the number of all alleles found in list: the number of the subgroup/gene/allele is divided by the number of all alleles mentioned for all sequences (in the case there are more alleles/genes mentioned for one sequence).

The \code{PDF} character string should be only the project name (without ".pdf"). For an analysis of the gene usage distribution, a plot called "PDF"_Gene-usage.pdf will be saved to the working direktory. 

If \code{geneUsage.Func = T} gene usage vs. functionality will be analyzed (and a figure called "PDF"_Gene-usage_vs_Functionality.pdf will be saved) and if \code{geneUsage.JunctionFr = T} gene usage vs. junction frame usage will be analyzed (and figure called "PDF"_Gene-usage_vs_Junction-frame.pdf will be saved).
}
\value{
Output is a list containing
     \item{gene_usage}{relative or absolute abundances of subgroups/genes}
     \item{gene_usage_vs_functionality}{(optional)}
     \item{gene_usage_vs_junction_frame}{(optional)}
}

\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
\code{\link{geneUsage}}, \code{\link{plotGeneUsage}}
}
\examples{
data(clones.ind)
geneUsage(genes = clones.ind$V_gene, level = "subgroup", 
     functionality = clones.ind$Functionality_all_sequences)
\dontrun{plotGeneUsage(geneUsage.tab = clones.ind$V_gene, 
     plotFunctionality = TRUE)}
}

