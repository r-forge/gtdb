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
\name{qqprcomp}
\alias{qqprcomp}
\title{Quantile-Quantile Plots of PCA Loadings}
\description{
  Generates a series of quantile-quantile plots of either sample or
  SNP loadings from a principal components analysis.
}
\usage{
qqprcomp(x, col=1:6, layout, ...)
}
\arguments{
  \item{x}{either a PCA result structure returned by
    \code{\link{prcomp.gt.data}} or \code{\link{prcomp.gt.dataset}},
    or a table of SNP loadings from \code{\link{snp.loadings}}.}
  \item{col}{a vector of component numbers to plot.}
  \item{layout}{a vector with two elements giving the numbers of
    columns and rows to use for arranging the plot panels.  If not
    specified, the plots are arranged to minimize wasted space.}
  \item{\dots}{additional arguments passed to \code{\link{qqthin}}.}
}
\details{
  This function generates a series of normal Q-Q plots of loadings
  from a principal components analysis.  It accepts either a PCA
  result structure (for sample loadings), or a table of SNP loadings
  returned by \code{\link{snp.loadings}}.  Sparse Q-Q plots are
  generated automatically when the number of loadings is very large.

  Q-Q plots are useful for distinguishing components that reveal
  either sample or genomic structure in the data, from components
  that do not.  If a component has flat Q-Q plots for both sample and
  SNP loadings, it is unlikely to be informative.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{prcomp.gt.data}}, \code{\link{prcomp.gt.dataset}},
  \code{\link{snp.loadings}}, \code{\link{qqthin}}.
}
\examples{
gt.demo.check()
p1 <- prcomp(gt.dataset('Demo_1', gt.filter=(ploidy=='A')))
qqprcomp(p1)
pt <- fetch.pt.data('Demo_1')
p2 <- prcomp(gt.dataset('Demo_1', gt.filter=(ploidy=='A')),
             sample.mask=is.na(pt$father))
qqprcomp(p2)
}
\keyword{hplot}
