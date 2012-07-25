% $Id:$
\name{exactTestafterTbT}
\alias{exactTestafterTbT}
\title{exactTest after TMM-baySeq-TMM procedure}
\description{
  This function perform exact test with edgeR after TMM-baySeq-TMM procedure
  via \code{\link{do_TbT}}.
}
\usage{
   exactTestafterTbT(names, counts, group, sample_num=10000, prop.used=NULL, span=NULL, grid.length=500)
}
\arguments{
\item{names}{A vector containing the name of each element eg gene}
\item{counts}{The data matrix to be analysed. Numerical data only}
\item{group}{A vector describing the data class for columns in data}
\item{sample_num}{Sample number for Bayesian estimation}
\item{prop.used}{parameter passed to estimate dispersion (default 0.5 for edgeR =< 2.6)}
\item{span}{parameter passed to estimate dispersion (changed from prop.used after edgeR 2.7)}
\item{grid.length}{parameter passed to estimate dispersion}
}
\value{
    This function return a list of data, norm_f_TbT, counts, and group. The data is a 
    data frame containing the names and counts
    as the input, table output by the exactTest(), FDR calculated by
    p.adjust(), and rank by the FDR. norm_f_TbT is calculated by do_TbT().
    The counts and group are the copy of the input.

    \item{names}{As in the input argument.}
    \item{counts}{As in the input argument.}
    \item{table}{As returned from exactTest().}
    \item{FDR}{False discovery rate calculated by p.adjust().}
    \item{rank_edgeR}{The rank of the above.}
}
\seealso{
\code{\link{do_TbT}}
\code{\link{exactTest}}
\code{\link{p.adjust}}
}

\examples{
  sample <- NBsample()
  \dontrun{
  out <- exactTestafterTbT(paste("gene",1:nrow(sample$data), sep=""), sample$data, c(1,1,1,2,2,2)) }
  \dontshow{out <- exactTestafterTbT(paste("gene",1:nrow(sample$data), sep=""), sample$data, c(1,1,1,2,2,2), 10)}
}
% vim: syntax=tex
