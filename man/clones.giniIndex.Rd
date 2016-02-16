\name{clones.giniIndex}
\alias{clones.giniIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gini index of clones
}
\description{
This function calculates the Gini index of clones.
}
\usage{
clones.giniIndex(clone.size = NULL, PDF = NULL)
}

\arguments{
  \item{clone.size}{
     Vector containing clone sizes
}
  \item{PDF}{
     PDF project name (see Details)
}
}
\details{
This function calculates the Gini index of clones. Input is a vector containing clone sizes (copy number). The Gini index measures the inequality of clone size distribution. It's between 0 and 1. An index of 0 represents a polyclonal distribution, where all clones have same size. An index of 1 represents a perfect monoclonal distribution.

The \code{PDF} character string should be only the project name (without ".pdf").    
A figure called "PDF"_Lorenz-curve.pdf will be saved to the working directory. The Lorenz curve represents the clone distribution.

}

\references{
Cowell, F. A. (2000): Measurement of Inequality. in A B Atkinson/F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam

Gastwirth, J. L. (1972): The Estimation of the Lorenz Curve and Gini Index. The Review of Economics and Statistics 54 (3): 306-316

Lorenz, M. O. (1905): Methods of measuring the concentration of wealth. Publications of the American Statistical Association 9 (70): 209-219

Zeileis, A. (2014): ineq: Measuring Inequality, Concentration, and Poverty. R package version 0.2-13. http://CRAN.R-project.org/package=ineq
}
\author{
Julia Bischof
}

\examples{
data(clones.ind)

gini<-clones.giniIndex(clone.size=clones.ind$total_number_of_sequences, PDF = "Example")
}
