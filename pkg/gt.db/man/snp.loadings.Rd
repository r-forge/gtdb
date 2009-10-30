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
\name{snp.loadings}
\alias{snp.loadings}
\title{SNP Loadings from Principal Components Analyses}
\description{
  Computes loadings of assays onto previously computed principal
  components.  Results can be computed for either a specific set of
  SNPs or the entire original genotype dataset.
}
\usage{
snp.loadings(x, data)
}
\arguments{
  \item{x}{a structure returned by \code{\link{prcomp}}.}
  \item{data}{either a data frame of genotypes from
    \code{\link{fetch.gt.data}}, or a dataset description from
    \code{\link{gt.dataset}}.} 
}
\details{
  \code{snp.loadings} computes SNP loadings across either an entire
  dataset or a specific set of genotype assays.  The genotype data
  needs to cover the same individuals included in the original
  principal components analysis.

  If the analysis was done with \code{\link{prcomp.gt.dataset}} and
  \code{data} is missing, then loadings are computed across the same
  dataset specification used in the original analysis.
}
\value{
  A list with the following components:
  \item{loadings}{a data frame with one row per assay, and \code{nc}
    columns, containing the loadings of each assay onto the top
    principal components.}
  \item{sdev}{the standard deviations of the top \code{nc} principal
    components.}
  \item{assays}{annotations for the assays in \code{loadings}.}
  \item{dataset.name}{the name of the source genotype dataset.}
  \item{call}{the call used to generate the analysis.}
}
\seealso{
  \code{\link{prcomp.gt.data}},
  \code{\link{fetch.gt.data}},
  \code{\link{gt.dataset}},
  \code{\link{apply.loadings}}.
}
\examples{
gt.demo.check()
p <- prcomp(gt.dataset('Demo_1', gt.filter=(ploidy=='A')))
s <- snp.loadings(p, gt.dataset('Demo_1'))
p$loadings[1:5,1:5]
s$loadings[1:5,1:5]
}
\keyword{multivariate}
