\name{sequences.distance}
\alias{sequences.distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity/distance indices for sequence data
}
\description{
This function calculates different dissmilarity/distance indices of sequences.
}
\usage{
sequences.distance(sequences = NULL, groups = NULL, 
     method = c("levenshtein", "cosine", "q-gram", "jaccard", "ja-wi", 
                "dam-le", "hamming", "osa", "lcs"), divLength = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{
     Vector containing sequences
}
  \item{groups}{
     Vector containing names of different samples (if present)
}
  \item{method}{
     Dissmilariy method (see details)
}
  \item{divLength}{
     Divide sequences into subsets of the same sequence length? (default: FALSE)
}
}
\details{
This function calculates dissmiliarity/distance indices based on sequences. Levenshtein, cosine, q-gram, Jaccard, Jaro-Winker (\code{ja-wi}), Damerau-Levenshtein (\code{dam-le}), Hamming, Optimal string alignment (\code{osa}) and longest common substring (\code{lcs}) distance can be chosen. For details see \code{\link[stringdist]{stringdist-metrics}}.
}
\value{
Output is a distance matrix containing dissimilarity indices/distances between sequences.
}
\references{
van der Loo M (2014). The stringdist package for approximate string matching. The R Journal, 6, pp. 111-122. \url{https://CRAN.R-project.org/package=stringdist}
}
\author{
Julia Bischof
}

\seealso{
\code{\link{dist.PCoA}}, \code{\link{plotDistPCoA}}, \code{\link{geneUsage.distance}}
}
\examples{
\dontrun{
data(clones.ind)
data(clones.allind)
dist1<-sequences.distance(sequences = clones.ind$unique_CDR3_sequences_AA, 
     method = "levenshtein", divLength=TRUE)
dist2<-sequences.distance(sequences = clones.allind$unique_CDR3_sequences_AA, 
     groups = clones.allind$individuals, method = "cosine", divLength=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
