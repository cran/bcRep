\name{readIMGT}
\alias{readIMGT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read IMGT/HighV-QUEST output files
}
\description{
This function reads IMGT/HighV-QUEST output files into a data frame.
}
\usage{
readIMGT(data, filterNoResults = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     path and file name of IMGT/HighV-QUEST output file
}
  \item{filterNoResults}{
     Shall sequences without any information be excluded? (default: TRUE)
}
}
\details{
"-" or spaces are be converted to "_" (single quotes are ignored).
}
\value{
Output is a data frame containing IMGT/HighV-QUEST output data.
}
\references{
IMGT/V-QUEST Documentation: \url{http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#output3}
}
\author{
Julia Bischof
}

\examples{
\dontrun{
tab<-readIMGT("PathTo/file.txt",filterNoResults=TRUE)
str(tab)
}
}

