\name{simulateReadCounts}
\alias{simulateReadCounts}
\title{Generate simulation data from negative binomial (NB) distribution}
\description{
This function generates simulation data with a specified condition. 
It can generate not only all of the simulation data analyzed in Kadota 
et al., (2012) but also simulation data with more complex design 
such as two or more groups and/or without replicates.
}
\usage{
simulateReadCounts(Ngene = 10000, PDEG = 0.20, DEG.assign = c(0.9, 0.1),
                  DEG.model = NULL, DEG.foldchange = NULL,
                  replicates = c(3,3))
}
\arguments{
	\item{Ngene}{numeric scalar specifying the number of genes.}
	\item{PDEG}{numeric scalar specifying the proportion of
	 differentially expressed genes (DEGs).}
	\item{DEG.assign}{numeric vector specifying the proportions
	 of DEGs up-regulated in individual groups to be compared.
	 The number of element should be the same as the number of
	 groups compared.}
	\item{DEG.model}{character string specifying the distribution
	 for sampling the fold-change (FC) values. Possible values 
 	 are \code{"uniform"} (default) and \code{"gamma"}.}
	\item{DEG.foldchange}{list. If \code{DEG.model = "uniform"},
	 \eqn{i}-th component contains numeric scalar specifying the
	 degree of FC for Group \eqn{i}. The default is
	 \code{DEG.foldchange = list(4, 4)}, indicating that the levels
	 of DE are four-fold in both groups. If \code{DEG.model = "gamma"},
	 \eqn{i}-th component contains numeric vector consisting of three
	 elements, i.e., c(\eqn{x}, \eqn{y}, \eqn{z}), for determining
	 the FC value for DEGs up-regulated in Group \eqn{i}. Specifically,
	 the FC value for each DEG is randomly sampled from
	 "\eqn{x} + a gamma distribution with shape (= \eqn{y}) and
	 scale (= \eqn{z}) parameters", indicating that the minimum and
	 average FC values for DEGs correspond to \eqn{x} and
	 (\eqn{x} + \eqn{y} * \eqn{z}).}
	\item{replicates}{numeric vector indicating the numbers of (biological) 
	replicates for individual groups compared.}
}
\details{
The empirical distribution of read counts 
used in this function is calculated from a RNA-seq dataset 
obtained from \emph{Arabidopsis} data 
(three biological replicates for both the treated and non-treated samples), 
the arab object, in NBPSeq package (Di et al., 2011). 
The overall design about the simulation conditions introduced 
can be viewed as a pseudo-color image by the \code{\link{plotFCPseudocolor}} function.
}
\value{
A \link{TCC-class} object containing following fields:
	\item{count}{numeric matrix of simulated count data.}
	\item{group}{numeric vector indicating the experimental
	 group for each sample.}
	\item{norm.factors}{numeric vector as a placeholder for
	 normalization factors.}
	\item{stat}{list for storing results after the execution of
	 the \code{\link{calcNormFactors}} (and \code{\link{estimateDE}})
	 function.}
	\item{estimatedDEG}{numeric vector as a placeholder for indicating
	 which genes are up-regulated in particular group compared to the
	 others. The values in this field will be populated after the
	 execution of the \code{\link{estimateDE}} function.}
	\item{simulation}{list containing four fields: \code{trueDEG},
	 \code{DEG.foldchange}, \code{PDEG}, and \code{group}. The
	 \code{trueDEG} field (numeric vector) stores information about
	 DEGs: 0 for non-DEG, 1 for DEG up-regulated in Group 1, 2 for
	 DEG up-regulated in Group 2, and so on. The information for
	 the remaining three fields is the same as those indicated in
	 the corresponding arguments.}
}
\examples{
# Generating a simulation data for comparing three groups
# (G1 vs. G2 vs. G3) with biological replicates. 
# the first 30% genes are DEGs, where the 70%, 20%, and 10% are
# up-regulated in G1, G2, G3, respectively. The levels of DE are
# 3-, 10, and 6-fold in individual groups.
tcc <- simulateReadCounts(Ngene = 100, PDEG = 0.3, 
                         DEG.assign = c(0.7, 0.2, 0.1),
                         DEG.foldchange = c(3, 10, 6), 
                         replicates = c(2, 4, 3))
dim(tcc$count)
head(tcc$count)
str(tcc$simulation)
head(tcc$simulation$trueDEG)
}
\keyword{methods}
