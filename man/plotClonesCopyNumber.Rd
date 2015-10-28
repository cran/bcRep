\name{plotClonesCopyNumber}
\alias{plotClonesCopyNumber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copy number of clones
}
\description{
This Function plots the copy number distribution of clones.
}
\usage{
plotClonesCopyNumber(copyNumber = NULL, withOutliers = TRUE, color = "gray", 
     title = NULL, PDF=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copyNumber}{
     Vector of copy numbers of clones
}
  \item{withOutliers}{
     Shall outliers (anything bigger than 75% quantile) be included? (default: TRUE)
}
\item{color}{
	color used for plots (default: gray)
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
The \code{PDF} character string should be only the project name (without ".pdf"). 

A figure called "PDF"_Clone_copy_number.pdf will be saved to your working directory.
}
\author{
Julia Bischof
}
\seealso{
 \code{\link{clones.CDR3Length}}, \code{\link{plotClonesCDR3Length}}
}
\examples{
data(clones.ind)
\dontrun{plotClonesCopyNumber(copyNumber = clones.ind$total_number_of_sequences, 
     color = "red", title = "Copy number distribution")}
}

