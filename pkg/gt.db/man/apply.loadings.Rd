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
\name{apply.loadings}
\alias{apply.loadings}
\title{Apply PCA SNP Loadings to a New Genotype Set}
\description{
  Takes a table of SNP loadings from a principal components analysis,
  and applies them to another set of genotypes, yielding a new set of
  sample loadings.
}
\usage{
apply.loadings(x, data)
}
\arguments{
  \item{x}{SNP loadings from \code{\link{snp.loadings}}.}
  \item{data}{either a data frame of genotype data from
    \code{\link{fetch.gt.data}}, or a dataset description from
    \code{\link{gt.dataset}}.} 
}
\details{
  If the SNP loadings were computed on the same platform as the new
  dataset, then corresponding assays are identified by name.  Across
  platforms, assays are matched up by genomic position.
}
\value{
  A data frame of sample loadings, with one column per principal
  component, with samples ordered as they are in the supplied genotype
  strings.
}
\seealso{
  \code{\link{prcomp}},
  \code{\link{snp.loadings}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_1')
# compute PCA for unrelated individuals
p1 <- prcomp(gt.dataset('Demo_1'), is.na(pt$father))
# now compute sample loadings for all individuals
s <- snp.loadings(p1)
p2 <- apply.loadings(s, gt.dataset('Demo_1'))
}
\keyword{manip}
