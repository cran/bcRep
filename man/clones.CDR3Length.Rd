\name{clones.CDR3Length}
\alias{clones.CDR3Length}
\alias{plotClonesCDR3Length}
\title{
CDR3 length distribution of clones
}
\description{
This function gives information about the CDR3 length distribution of clones. Results can be returned as relative or absolute values and be visualized as a barplot, using \code{plotClonesCDR3Length}.
}
\usage{
clones.CDR3Length(CDR3Length = NULL, functionality = NULL, junctionFr = NULL, 
     abundance=c("relative","absolute"), ...)
     
plotClonesCDR3Length(CDR3Length=NULL,functionality=NULL, junctionFr=NULL,   
     color=c("orange","darkblue","gray"), abundance=c("relative","absolute"), 
     title=NULL, PDF=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CDR3Length}{
      Vector of CDR3 length of clones (amino acids or nucleotides)
}
  \item{functionality}{
     Vector of functionality of clones (same order than CDR3Length)
}
  \item{junctionFr}{
     Vector of junction frame usage of clones (same order than CDR3Length)
}
  \item{abundance}{
     Shall relative or absolute values be returned? (default: relative)
}
\item{color}{
	color used for plots (default: c("orange","darkblue","gray"))
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
The \code{PDF} character string should be only the project name (without ".pdf"). For simple analysis of the CDR3 length, a plot called "PDF"_CDR3-length.pdf will be saved to the working direktory. 

If \code{CDR3Length.Func = T} a figure containing CDR3 length vs. functionality is saved (called "PDF"_CDR3-length_vs_Functionality.pdf) and 
if \code{CDR3Length.JunctionFr = T} a figure containing CDR3 length vs. junction fram usage (called "PDF"_CDR3-length_vs_Junction-frame.pdf) is saved.
}
\value{
Output is a list containing
     \item{CDR3_length}{relative or absolute abundances of CDR3 length}
     \item{CDR3_length_vs_functionality}{(optional)}
     \item{CDR3_length_vs_junction_frame}{(optional)}
}

\author{
Julia Bischof
}

\seealso{
 \code{\link{clones.CDR3Length}}, \code{\link{plotClonesCDR3Length}}, \code{\link{plotClonesCopyNumber}}, \code{\link{geneUsage}}
}
\examples{
\dontrun{data(clones.ind)
clones.CDR3Length(CDR3Length = clones.ind$CDR3_length_AA, 
     functionality = clones.ind$Functionality_all_sequences)
plotClonesCDR3Length(CDR3Length = clones.ind$CDR3_length_AA, 
     functionality = clones.ind$Functionality_all_sequences)}
}

