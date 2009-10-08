%
% Copyright (C) 2009, Perlegen Sciences, Inc.
% 
% Written by David A. Hinds <dhinds@sonic.net>
% 
% This is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the license, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>
% 
\name{prcomp.gt.data}
\alias{prcomp.gt.data}
\alias{prcomp.gt.dataset}
\title{Principal Components Analysis of Genotype Data}
\description{
  Performs a principal components analysis of a matrix of genotypes
  and scores each sample against the most significant components.
}
\usage{
\method{prcomp}{gt.data}(x, sample.mask=TRUE, nc=20, ...)
\method{prcomp}{gt.dataset}(x, sample.mask=TRUE, nc=20, ...)
}
\arguments{
  \item{x}{a data frame of genotypes from \code{\link{fetch.gt.data}},
    or a dataset description from \code{\link{gt.dataset}}.}
  \item{sample.mask}{a logical vector for subsetting samples.}
  \item{nc}{the number of components to return.}
  \item{\dots}{not used.}
}
\details{
  This function is useful for computing principal components in cases
  where an entire genotype matrix can be accommodated in memory.
  Genotype vectors are centered and scaled, and missing values are
  imputed, as in Price \emph{et al.} (2006).  The principal components
  are computed using \code{\link[stats:prcomp]{prcomp}}.
}
\value{
  A list with the following components:
  \item{loadings}{a data frame with one row per row in the \code{data}
  table, and \code{nc} columns, containing the loadings of each sample
  onto the top principal components.}
  \item{sdev}{the standard deviations of the top \code{nc} principal
  components.}
  \item{assays}{the number of assays included in the analysis.}
  \item{dataset.name}{the name of the genotype dataset.}
  \item{call}{the call used to generate the analysis.}
}
\references{
  Price, A. L., \emph{et al.} (2006) Principal components analysis
  corrects for stratification in genome-wide association studies.
  \emph{Nat. Genet.} \bold{38}: 904-909.
}
\seealso{
  \code{\link[stats:prcomp]{prcomp}},
  \code{\link{fetch.gt.data}},
  \code{\link{gt.dataset}},
  \code{\link{snp.loadings}},
  \code{\link{apply.loadings}}.
}
\examples{
gt.demo.check()
g <- fetch.gt.data('Demo_1')
p <- prcomp(subset(g, (ploidy=='A')))
screeplot(p)
with(p$loadings, plot(PC1, PC2))
}
\keyword{multivariate}
