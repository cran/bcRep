\name{compare.aaDistribution}
\alias{compare.aaDistribution}
\alias{plotCompareAADistribution}
\title{
Compare amino acid distribution of different samples
}
\description{
This function compares the amino acid distribution of different samples. Sequences of the same length are clustered and analyzed. Additionally the number of sequences for each sample can be returned. Parallel processing is possible.
}
\usage{
compare.aaDistribution(sequence.list = NULL, names = NULL, numberSeq = FALSE, 
     nrCores = 1)

plotCompareAADistribution(comp.tab = NULL, plotSeqN = FALSE, colors = NULL, 
     title = NULL, PDF = NULL)
}

\arguments{
  \item{sequence.list}{
     A list containing vectors of amino acid sequences of each sample
}
  \item{names}{
     A vector containing names for the samples (default: Sample1, Sample2, ...)
}
  \item{numberSeq}{
     Shall number of sequences used for analysis be returned? (default: FALSE)
}
  \item{nrCores}{
     Number of cores used for parallel processing
}
  \item{comp.tab}{
     Output tab from \code{compare.aaDistribution()}
}
  \item{plotSeqN}{
     Shall number of sequences used for analysis be plotted? (only possible, if in \code{compare.aaDistribution(...,numberSeq=TRUE,...)} is used; default: FALSE)
}
  \item{colors}{
     colors used for individuals in sequence number plot (default: rainbow)
}
  \item{title}{
     Title of plot
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
Amino acid distribution for each individual is analyzed. Sequenzes of the same length are clustered and amino acid distribution for each position is measured.

The \code{PDF} character string should be only the project name (without ".pdf"). 

A figure called "PDF"_Comparison_Amino-acid-distribution.pdf will be saved to the working directory. 

Optional another figure called PDF"_Comparison_Number-of-sequences.pdf" will be returned.

}
\value{
Output is a list containing amino acid distributions for each sequence length and each sample.
}
\author{
Julia Bischof
}
\note{
     For large datasets computational time can be extensive.
}

\seealso{
\code{\link{compare.aaDistribution}}, \code{\link{plotCompareAADistribution}}, \code{\link{aaDistribution}}
}
\examples{
data(aaseqtab)
data(aaseqtab2)

AAdistr.comp<-compare.aaDistribution(sequence.list = list(aaseqtab$CDR3_IMGT, 
     aaseqtab2$CDR3_IMGT), names = c("IndA", "IndB"), numberSeq = FALSE, nrCores = 1)
\dontrun{
plotCompareAADistribution(comp.tab = AAdistr.comp, plotSeqN = FALSE, PDF = "Example")
}
}

