\name{geneUsage.distance}
\alias{geneUsage.distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity/distance indices for gene usage data
}
\description{
This function calculates Bray-Curtis, Jaccard or cosine indices for gene usage data of different samples.
}
\usage{
geneUsage.distance(geneUsage.tab=NULL, names=NULL, 
     method=c("bc","jaccard", "cosine"), cutoff=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneUsage.tab}{
     gene usage table with genes as columns and samples as rows
}
  \item{names}{
     Names of samples (default: Sample1...n)
}

  \item{method}{
     Distance/dissimilarity index to be used for calculation. On of Bray-Curtis (bc), Jaccard (jaccard) or cosine (cosine)
}
  \item{cutoff}{
     Cutoff for gene proportions (default: 0; see details)
}
}
\details{
This function calculates dissmiliarity indices based on gene usage data of different samples (columns = genes, rows = samples). Bray-Curtis, Jaccard or cosine indices can be chosen.

For explanation of Bray-Curtis and Jaccard index see \code{\link[vegan]{vegdist}}.

For explanation of cosine index see \code{\link[proxy]{dist}}.

When using Jaccard index, a cutoff for gene proportions can be given. Proportions will be transformed into absence/presence data (<= cutoff; > cutoff) and afterwards Jaccard indices are calculated.
}
\value{
Output is a matrix containing dissimilarity/distance indices between samples, based on gene usages.
}
\references{
Bray, J. R. and J. T. Curtis. 1957. An ordination of upland forest communities of southern Wisconsin. Ecological Monographs 27:325-349.

Levandowsky, Michael; Winter, David (1971), "Distance between sets", Nature 234 (5): 34-35

Graham L. Giller (2012). "The Statistical Properties of Random Bitstreams and the Sampling Distribution of Cosine Similarity". Giller Investments Research Notes (20121024/1)

Jari Oksanen, F. Guillaume Blanchet, et al. (2015). vegan: Community Ecology Package. R package version 2.3-1. \url{https://CRAN.R-project.org/package=vegan}
  
David Meyer and Christian Buchta (2015). proxy: Distance and Similarity Measures. R package version 0.4-15. \url{https://CRAN.R-project.org/package=proxy}
}
\author{
Julia Bischof
}
\seealso{
 \code{\link{dist.PCoA}},  \code{\link{plotDistPCoA}}, \code{\link{sequences.distance}}
}
\examples{
data(vgenes) # VH gene proportions of 10 samples (rows) and 30 VH genes (columns)
disttab<-geneUsage.distance(geneUsage.tab = vgenes, method = "bc")

}

