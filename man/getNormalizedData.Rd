\name{getNormalizedData}
\alias{getNormalizedData}
\title{Obtain normalized count data}
\description{
This function generates normalized count data from both original count
data and calculated normalization factors. 
}
\usage{getNormalizedData(tcc)}
\arguments{
	\item{tcc}{\link{TCC-class} object.}
}
\details{
This function is generally used after the \code{\link{calcNormFactors}}
function that calculates normalization factors. 
The normalized data is calculated using both the original count data 
stored in the \code{count} field and the normalization factors
stored in the \code{norm.factors} field in the \link{TCC-class} object.
}
\value{
A numeric matrix containing normalized count data.
}
\examples{
\dontrun{
# Note that the hypoData has non-DEGs at 2001-10000th rows
nonDEG <- 2001:10000
data(hypoData)
summary(hypoData[nonDEG, ])

# Obtaining normalized count data after performing the 
# TbT normalization method (Kadota et al., 2012) from hypoData,
# i.e., TbT-normalized data
group <- c(3, 3)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq")
normalized.count <- getNormalizedData(tcc)
summary(normalized.count[nonDEG, ])

# Obtaining normalized count data after performing the 
# DEGES/edgeR normalization method, i.e., DEGES/edgeR-normalized data
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc)
normalized.count <- getNormalizedData(tcc)
summary(normalized.count[nonDEG, ])

# Obtaining normalized count data after performing the 
# TMM normalization method (Robinson and Oshlack, 2010), 
# i.e., TMM-normalized data
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", iteration = FALSE)
normalized.count <- getNormalizedData(tcc)
summary(normalized.count[nonDEG, ])
}
}
\keyword{methods}
