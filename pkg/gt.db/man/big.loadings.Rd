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
\name{big.loadings}
\alias{big.loadings}
\title{Identify Large Sample or SNP Loadings from a PCA Analysis}
\description{
  Takes either a PCA result, or SNP loadings from PCA, and returns a
  data frame describing large loadings for each principal component.
}
\usage{
big.loadings(x, sigma=6)
}
\arguments{
  \item{x}{either a PCA result structure, or SNP loadings.}
  \item{sigma}{number of standard deviations to qualify as \dQuote{big}.}
}
\value{
  A data frame with sample or assay identifiers, and the following
  additional columns:
  \item{name}{the principal component with high loading.}
  \item{value}{the loading for this sample and component.}
  \item{zscore}{the standardized loading for this sample.}
  \item{var}{the fraction of variance in this component accounted
    for by this sample.}
}
\seealso{
  \code{\link{prcomp.gt.data}},
  \code{\link{prcomp.gt.dataset}},
  \code{\link{snp.loadings}}.
}
\examples{
gt.demo.check()
g <- fetch.gt.data('Demo_1')
pc <- prcomp(subset(g, (ploidy=='A')))
big.loadings(pc,5)
sl <- snp.loadings(pc, g)
big.loadings(sl,5)
}
\keyword{manip}
