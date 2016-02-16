\name{dist.PCoA}
\alias{dist.PCoA}
\alias{plotDistPCoA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal coordinate analysis (PCoA; multidimensional scaling [MDS]) of dissimilarity/distance indices
}
\description{
Performs and plots a principal coordinate analysis (PCoA) of dissimilarity/distance indices. Correction methods can be used. Merging of samples to groups is possible in the plot function.
}
\usage{
dist.PCoA(dist.tab = NULL, correction = c("lingoes", "cailliez", "none"))

plotDistPCoA(pcoa.tab = NULL, groups = NULL, names = NULL, axes = NULL, 
     plotCorrection = FALSE, title = NULL, plotLegend=FALSE, PDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.tab}{
     Dissimilarity/distance matrix (e.g. from \code{sequences.distance()})
}
  \item{correction}{
     Correction method of PCoA: Lingoes, Cailliez or none
}

  \item{pcoa.tab}{
     PCoA outout from \code{dist.PCoA()}
}
  \item{groups}{
     data.frame containing sequences (1. column) and groups (2. column)
}
  \item{names}{
     Names of samples/axes
}
  \item{axes}{
     Which axes shall be plotted? e.g. c(1,2) for axes 1 and 2
}
  \item{plotCorrection}{
     Shall corrected or uncorrected eigenvalues be plotted?
}
  \item{title}{
     Title of the plot
}
  \item{plotLegend}{
     Shall legend be plotted?
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
This function provides a PCoA object for dissimilarity indices/distances as input (e.g. from functions \code{sequences.distance} or \code{geneUsage.distance()}). For further details of pcoa see \code{\link[ape]{pcoa}}.

The plot function provides a figure with the principal coordinates with positive eigenvalues (in the case of no correction) or the principal coordinates with positive eigenvalues from the distance matrix corrected using the specified correction method. The principal coordinates correspond to the specified axes.

A figure called "PDF"_PCoA.pdf will be saved to the working directory.
}
\value{
Output is an PCoA object, see \code{\link[ape]{pcoa}}.

}
\references{
Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
}
\author{
Julia Bischof
}

\seealso{
\code{\link{dist.PCoA}}, \code{\link{plotDistPCoA}}, \code{\link{sequences.distance}},    
\code{\link{geneUsage.distance}}, \code{\link[ape]{pcoa}}
}
\examples{
\dontrun{
data(clones.ind)
seq.dist<-sequences.distance(sequences = clones.ind$unique_CDR3_sequences_AA, 
     method = "levenshtein", divLength=F)
distpcoa<-dist.PCoA(dist.tab = seq.dist, correction = "none")

# 'groups' data.frame for plot function: in the case, there are no groups:
groups.vec<-unlist(apply(data.frame(clones.ind$unique_CDR3_sequences_AA),1,
          function(x){strsplit(x,split=", ")[[1]]}))
groups.vec<-cbind(groups.vec, 1)

plotDistPCoA(pcoa.tab = distpcoa, groups=groups.vec, axes = c(1,2), 
     plotCorrection = FALSE, title = NULL, plotLegend=T, PDF = "TEST")    
}
}
