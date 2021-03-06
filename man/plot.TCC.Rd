\name{plot}
\docType{methods}
\alias{plot.TCC}
\alias{plot}
\title{Plot an log fold-change versus log average expression
 (so-called M-A plot)}
\description{This function generates a scatter plot of log fold-change
 (i.e., \eqn{M = \log_2 G2 - \log_2 G1}{M = log2(G2) - log2(G1)} on 
 the \eqn{y}-axis) between Groups 1 vs. 2) versus log average 
 expression (i.e., \eqn{A = (\log_2 G1 + \log_2 G2) / 2}{A = (log2(G1) + 
 log2(G2)) / 2} on the \eqn{x}-axis) using normalized count data.}
\usage{
\method{plot}{TCC}(x, FDR = NULL, median.lines = FALSE, floor = 0,
     main = NULL,
     xlab = expression(A == (log[2] * G2 + log[2] * G1) / 2),
     ylab = expression(M == log[2] * G2 - log[2] * G1),
     xlim = NULL, ylim = NULL, cex = 0.3, pch = 19, col = NULL, ...)
}
\arguments{
	\item{x}{\link{TCC-class} object.}
	\item{FDR}{numeric scalar specifying a false discovery rate (FDR)
	 threshold for determining differentially expressed genes (DEGs)}
	\item{median.lines}{logical. If \code{TRUE}, horizontal lines
	 specifying the median M values for non-DEGs (black), DEGs
	 up-regulated in Group 1 (G1; blue), and G2 (red) are drawn.}
	\item{floor}{numeric scalar specifying a threshold for adjusting low count data.}
	\item{main}{character string indicating the plotting title.}
	\item{xlab}{character string indicating the \eqn{x}-label title.}
	\item{ylab}{character string indicating the \eqn{y}-label title.}
	\item{xlim}{numeric vector (consisting of two elements) specifying
	 the range of the \eqn{x} coordinates.}
	\item{ylim}{numeric vector (consisting of two elements) specifying
	 the range of the \eqn{y} coordinates.}
	\item{cex}{numeric scalar specifying the multiplying factor of the
	 size of the plotting points relative to the default (= 0.3).}
	\item{pch}{numeric scalar specifying a symbol or a single character
	 to be used as the default in plotting points.}
	\item{col}{vector specifying plotting color. The default is col 
	 = c(1, 4, 2, 5, 6, 7, ...) of color index.}
	\item{\dots}{further graphical arguments, see \code{\link{plot.default}}.}
}
\details{
This function generates roughly three different M-A plots 
depending on the conditions for \link{TCC-class} objects. 
When the function is performed just after the \code{new} method, 
all the genes (points) are treated as non-DEGs (the default is black; see Example 1). 
The \code{\link{simulateReadCounts}} function followed 
by the \code{\link{plot}} function can classify the genes as \emph{true} non-DEGs (black), 
\emph{true} DEGs up-regulated in Group 1 (G1; blue), 
and \emph{true} DEGs up-regulated in G2 (red) (see Example 2). 
The \code{\link{estimateDE}} function followed 
by the \code{\link{plot}} function generates \emph{estimated} DEGs (magenta) 
and the remaining \emph{estimated} non-DEGs (black).

Genes with normalized counts of 0 in any one group
cannot be plotted on the M-A plot because those M and A values 
cannot be calculated (as \eqn{\log 0} is undefined).  
Similar to the \code{\link[edgeR]{plotSmear}} function in edgeR package, 
\code{\link{plot}} function plots those points at the left side 
of the minimum A (i.e., log average expression) value. 
The \eqn{x} coordinate of those points is the minimum A value minus 1.
The \eqn{y} coordinate is calculated as if the 0 count was the minimum observed 
non-0 count in each group.
}
\value{
A scatter plot to the current graphic device. 
}
\examples{
# 1. 
# M-A plotting just after constructing the TCC class object from
# hypoData. In this case, the plot is generated from hypoData
# that has been scaled in such a way that the library sizes of 
# each sample are the same as the mean library size of the
# original hypoData. Note that all points are in black. This is
# because the information about DEG or non-DEG for each gene is 
# not indicated.
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
plot(tcc)

normalized.count <- getNormalizedData(tcc)
colSums(normalized.count)
colSums(hypoData)
mean(colSums(hypoData))

# 2-1. 
# M-A plotting of simulation data with default parameter settings.
# non-DEGs are in black, DEGs up-regulated in G1 and G2 are in
# blue and red.
\dontrun{
tcc <- simulateReadCounts()
plot(tcc, median.lines = TRUE)
}


# 2-2. 
# M-A plotting of DEGES/edgeR-normalized simulation data.
# It can be seen that the median M value for non-DEGs approaches
# zero.
\dontrun{
tcc <- simulateReadCounts()
tcc <- calcNormFactors(tcc)
plot(tcc, median.lines = TRUE)
}


# 3. 
# M-A plotting of DEGES/edgeR-normalized hypoData after performing
# DE analysis.
\dontrun{
data(hypoData)
group <- c(1, 1, 1, 2, 2, 2)
tcc <- new("TCC", hypoData, group)
tcc <- calcNormFactors(tcc)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
plot(tcc)
}


# Changing the FDR threshold
\dontrun{plot(tcc, FDR = 0.7)}
}
\keyword{methods}
