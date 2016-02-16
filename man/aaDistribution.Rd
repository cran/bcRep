\name{aaDistribution}
\alias{aaDistribution}
\alias{plotAADistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Amino acid distribution of sequences
}
\description{
This function calculates the amino acid distribution of sequences. Distribution is calculated for sequences of the same length and therein for each position. 

\code{aaDistribution} returns a list containing either only amino acid distribution or
a list containing amino acid distribution and analyzed number of sequences per length.

\code{plotAADistribution} visualizes the amino acid distribution of sequences of the same length.
}
\usage{
aaDistribution(sequences = NULL, numberSeq = FALSE)

plotAADistribution(aaDistribution.tab=NULL, plotSeqN=FALSE, 
     colors=NULL, PDF=NULL, ...)
}

\arguments{
  \item{sequences}{
	A vector containing amino acid sequences.
}
  \item{numberSeq}{
	TRUE: table containing number of sequences will be returned, as well (default: FALSE).
}
\item{aaDistribution.tab}{
	Output list of function \code{aaDistribution()}
}
\item{plotSeqN}{
	TRUE: Number of sequences for each length will be plotted (see Details; default: FALSE).
}
\item{colors}{
     Colors to be used for figure containing number of sequences (default: rainbow)
}
\item{PDF}{
     PDF project name (see Details)
}
\item{...}{

}
}
\details{
The vector containing sequences will be divided in sequences of the same length and then amino acid distribution for each position is analyzed.

If \code{numberSeq = T}, the number of sequences used for the analysis of sequences of the same length will be returned, as well. This information is also required for \code{plotAADistribution(...,    
plotSeqN = T)}. Sequence numbers equal to 0 are not plotted; the smallest number is 1.

The \code{PDF} character string should be only the project name (without ".pdf"). If \code{plotAADistr = T} a figure called "PDF"_Amino-acid-distribution.pdf will be saved to the working directory. If \code{plotSeqN = T} a figure called "PDF"_Number-of-sequences.pdf will be saved, as well.
}
\value{
Output is a list containing
\item{Amino_acid_distribution}{list contains data frames of amino acid distributions (including stop codons "*") for each length}
\item{Number_of_sequences_per_length}{data frame contains the number of sequences for each length, used for analysis (optional)}
}
\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}
\seealso{
\code{\link{trueDiversity}}
}
\examples{
data(aaseqtab)
aadistr<-aaDistribution(sequences = aaseqtab$CDR3_IMGT, numberSeq = TRUE)
\dontrun{plotAADistribution(aaDistribution.tab=aadistr, plotAADistr=TRUE, plotSeqN=FALSE, 
     PDF="test")}
}

