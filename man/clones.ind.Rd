\name{clones.ind}
\alias{clones.ind}
\docType{data}
\title{
Data frame: clones of one individual
}
\description{
This data frame represents the ouput of the function clones() for one individual (1000 clones).
}
\usage{data("clones.ind")}
\format{
  A data frame with 1000 observations on the following 14 variables.
  \describe{
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
data(clones.ind)
str(clones.ind)
boxplot(clones.ind$CDR3_length_AA,ylab="CDR3 length [AA]",xaxt="n",
     main="CDR3 length distribution of clones")
}
\keyword{datasets}
