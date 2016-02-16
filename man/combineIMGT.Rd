\name{combineIMGT}
\alias{combineIMGT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combination of several IMGT output folders
}
\description{
This function combines several IMGT ouput folders. IMGT/HighV-QUEST can analyse data with up to 500.000 sequences. In case one wants to analyse more than 500.000 sequences,  FASTA file have to be splitted into smaller files and have to be analysed individually. Afterwards IMGT output folders can be combined using this function.
}
\usage{
combineIMGT(folders = NULL, name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folders}{
     A list containing folder names
}
  \item{name}{
     Name of the new (combined) project
}
}

\value{
Output is a folder containing the 10 combined output files (no 11_Parameters.txt file).
}

\author{
Julia Bischof
}

\examples{
\dontrun{
## Combine folders IMGT1a, IMGT1b, IMGT1c to one folder named "NewProject"
combineIMGT(folders = c("pathTo/IMGT1a", "pathTo/IMGT1b", "pathTo/IMGT1c"), 
name = "NewProject)
}
}
