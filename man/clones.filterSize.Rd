\name{clones.filterSize}
\alias{clones.filterSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter clones for their size
}
\description{
This function filters clones for their size (clone copy number). It can be filtered for different kinds of tresholds (see Details).
}
\usage{
clones.filterSize(clones.tab = NULL, column = NULL, number = NULL, 
     propOfClones = NULL, propOfSequences = NULL, 
     method = c("two.tailed", "upper.tail", "lower.tail"))
}

\arguments{
  \item{clones.tab}{
     A data frame containing clones and their characteristics
}
  \item{column}{
     Name or number of column to be filtered
}
  \item{number}{
     An integer, giving treshold for size
}
  \item{propOfClones}{
     A proportion between 0 and 1, giving proportion of all clones (see Details)
}
  \item{propOfSequences}{
     A proportion between 0 and 1, giving proportion of all sequences (see Details)
}
  \item{method}{
     Filter method, either smallest and biggest numbers, only smallest or only biggest numbers (see Details)
}
}
\details{
This function filters clones for their size (clone copy number; total number of sequences belonging to a clone). It can be filtered for 

1) a given number of sequences, e.g. \code{number=20}: the 20 biggest and/or 20 smallest clones,

2) a proportion of all clones, e.g. \code{propOfClones=0.2}: the 20\% biggest and/or smallest clones,

3) a proportion of all sequences, e.g. \code{propOfSequences=0.01}: Clones, were more and/or less than 1\% of all sequences are included.

Only of these criteria is required.

\code{Method} determines which clones shall be returned: 1) biggest and smallest clones (\code{two.tailed}), 2) biggest clones (\code{upper.tail}) or 3) smallest clones (\code{lower.tail}). In case of \code{two.tailed} a list containing upper and lower tail will be returned.
}
\value{
Ouput is a data frame (\code{upper.tail, lower.tail}) or a list (\code{two.tailed}), containing filtered clones.
}

\author{
Julia Bischof
}
\seealso{
\code{\link{clones.filterFunctionality}}, \code{\link{clones.filterJunctionFrame}}
}
\examples{
data(clones.ind)

clones.filtered1<-clones.filterSize(clones.tab=clones.ind, 
     column="total_number_of_sequences", number=20, method="two.tailed")
clones.filtered2<-clones.filterSize(clones.tab=clones.ind, column=4, propOfClones=0.1, 
     method="upper.tail")
clones.filtered3<-clones.filterSize(clones.tab=clones.ind, column=4, 
     propOfSequences=0.02, method="lower.tail")
}
