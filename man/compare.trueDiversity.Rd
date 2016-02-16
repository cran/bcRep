\name{compare.trueDiversity}
\alias{compare.trueDiversity}
\alias{plotCompareTrueDiversity}
\title{
Compare richness or diversity of different samples
}
\description{
This function compares the richness or diversity of different samples. Analysis can be done for order 0 (richness), 1 (Shannon) and 2 (Simpson) (see Details). Parallel processing is possible.
}
\usage{
compare.trueDiversity(sequence.list = NULL, comp.aaDistribution.tab = NULL, 
     order = c(0, 1, 2), names = NULL, nrCores = 1)

plotCompareTrueDiversity(comp.tab = NULL, colors = NULL, title = NULL, PDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence.list}{
     A list containing vectors of amino acid sequences of each sample (see Details)
}
  \item{comp.aaDistribution.tab}{
     Output from \code{compare.aaDistribution()} (see Details)
}
  \item{order}{
     True diversity order (q). Values: 0, 1, 2
}
  \item{names}{
     A vector containing names for the samples (default: like Sample1, Sample2, ...)
}
  \item{nrCores}{
     Number of cores used for parallel processing
}
 \item{comp.tab}{
     Output tab from \code{compare.trueDiversity()}
}
  \item{colors}{
     Colors used for individuals (default: rainbow)
}
  \item{title}{
     Title of plot
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
This functions needs either a list containing vectors of sequences or the output of 

\code{compare.aaDistribution()} as input. In first case \code{compare.aaDistribution()} is first applied to data set. Richness or diversity is calculated for sequences of the same length, for each position. Analysis of true diversity of order \code{0, 1} and \code{2} is possible.    
Order 0: Richness (in this case it represents number of different amino acids per position).    
Order 1: Exponential function of Shannon entropy using the natural logarithm (weights all amino acids by their frequency).    
Order 2: Inverse Simpson entropy (weights all amino acids by their frequency, but weights are given more to abundant amino acids).    
These indices are very similar (Hill, 1973). For example the exponential function of Shannon index is linearly related to inverse Simpson.

\code{plotCompareTrueDiversity} returns an image with diversity plots for each length. Each plot contains the richness or diversity (y-axis) for each position (x-axis). Individuals are color coded.

The \code{PDF} character string should be only the project name (without ".pdf"). 

A figure called "PDF"_Comparison_True-diversity_q"order".pdf will be saved to the working directory. 
}
\value{
Output is a list containing 1) the diversity order, and 2) diversity values for each individual and each sequence length.
}
\references{
M. O. Hill: Diversity and Evenness: A Unifying Notation and Its Consequences; Ecology 54:2, p 427-432 (1973)

Lou Jost: Entropy and diversity; OIKOS 113:2 (2006)

Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens and Helene Wagner (2015). vegan: Community Ecology Package. R package version 2.3-0. http://CRAN.R-project.org/package=vegan
}
\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
\code{\link{compare.trueDiversity}}, \code{\link{plotCompareTrueDiversity}}, \code{\link{trueDiversity}}, 
\code{\link{compare.aaDistribution}}
}
\examples{
data(aaseqtab)
data(aaseqtab2)

trueDiv.comp<-compare.trueDiversity(sequence.list = list(aaseqtab$CDR3_IMGT, 
     aaseqtab2$CDR3_IMGT), names = c("IndA", "IndB"), order = 1, nrCores = 1)
\dontrun{
plotCompareTrueDiversity(comp.tab = trueDiv.comp, PDF = "Example")
}
}
