\name{clones.allind}
\alias{clones.allind}
\docType{data}
\title{
B cell clones of 8 individuals
}
\description{
This data frame represents the combined ouput of the function clones() for 7 individuals (in total 2948 clones).
}
\usage{data("clones.allind")}
\format{
  A data frame with 5286 observations on the following 15 variables.
  \describe{
    \item{\code{samples}}{a factor with 8 levels}
    \item{\code{unique_CDR3_sequences_AA}}{a character vector}
    \item{\code{CDR3_length_AA}}{a numeric vector}
    \item{\code{number_of_unique_sequences}}{a numeric vector}
    \item{\code{total_number_of_sequences}}{a numeric vector}
    \item{\code{sequence_count_per_CDR3}}{a character vector}
    \item{\code{V_gene}}{a character vector}
    \item{\code{V_gene_and_allele}}{a character vector}
    \item{\code{J_gene}}{a character vector}
    \item{\code{J_gene_and_allele}}{a character vector}
    \item{\code{D_gene}}{a character vector}
    \item{\code{all_CDR3_sequences_AA}}{a character vector}
    \item{\code{all_CDR3_sequences_nt}}{a character vector}
    \item{\code{Functionality_all_sequences}}{a character vector}
    \item{\code{Junction_frame_all_sequences}}{a character vector}
  }
}

\examples{
data(clones.allind)
str(clones.allind)
boxplot(clones.allind$CDR3_length_AA ~ clones.allind$samples,xlab="sample",
     ylab="CDR3 length [AA]",xaxt="n",main="CDR3 length distribution of clones", las=3)
axis(1,at=seq(1,length(unique(clones.allind$samples)),1),
     unique(clones.allind$samples))
}
\keyword{datasets}
