\name{trueDiversity}
\alias{trueDiversity}
\alias{plotTrueDiversity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
True diversity of sequences
}
\description{
This function provides information about the true diversity. Richness or diversity is calculated for sequences of the same length, for each position. Analysis of true diversity of order 0, 1 and 2 is possible. 

\code{trueDiversity} returns a list containing diversity indices.

\code{plotTrueDiversity} gives an overview about the richness or diversity of sequences with the same length.
}
\usage{
trueDiversity(sequences = NULL, aaDistribution.tab = NULL, order = c(0,1,2))

plotTrueDiversity(trueDiversity.tab=NULL,color="black", PDF=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{
	Vector containing sequences (see Details)
}
  \item{aaDistribution.tab}{
	Output of the function \code{AADistribution()} (see Details)
}
  \item{order}{
	True diversity order (q). Values: 0, 1, 2 (default: 1)
}
\item{trueDiversity.tab}{
	Output of function \code{trueDiversity()}
}
\item{color}{
	color used for plot (default: black)
}
\item{PDF}{
     PDF project name (see Details)
}
 \item{...}{

}
}
\details{
This functions needs either a vector of sequences or the output of \code{AADistribution()} as input. In first case \code{AADistribution()} is first applied to data set. Richness or diversity is calculated for sequences of the same length, for each position. Analysis of true diversity of order \code{0, 1} and \code{2} is possible. 
Order 0: Richness (in this case it represents number of different amino acids per position).
Order 1: Exponential function of Shannon entropy using the natural logarithm as the base (weights all amino acids by their frequency).
Order 2: Inverse Simpson entropy (weights all amino acids by their frequency, but weights are given more to abundant amino acids).
These indices are very similar (Hill, 1973). For example the exponential function of Shannon index is linearly related to inverse Simpson.

\code{plotTrueDiversity} returns an image with diversity plots for each length. Each plot contains the richness or diversity (y-axis) for each position (x-axis). 

The \code{PDF} character string should be only the project name (without ".pdf"). A figure called "PDF"_True-diversity_q"order".pdf will be saved to your working directory. 
}
\value{
Output is a list containing
  \item{True_diversity_order}{order of true diversity (q=0,1,2)}
  \item{True_diversity}{a list of true diversities for each position of each length}
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
     For large datasets computational time can be extensive for the calculation of amino acid proportions.
}
\seealso{
     \code{\link{aaDistribution}}, \code{\link{trueDiversity}}, \code{\link[vegan]{diversity}}
}
\examples{
data(aaseqtab)
trueDiv<-trueDiversity(sequences = aaseqtab$CDR3_IMGT, order = 1)
\dontrun{plotTrueDiversity(trueDiversity.tab=trueDiv,color="red")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

