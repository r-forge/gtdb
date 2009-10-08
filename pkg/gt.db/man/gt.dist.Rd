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
\name{gt.dist}
\alias{gt.dist}
\title{Calculate Pairwise Genotype Distances}
\description{
  Calculate genotype-based distances (broadly defined) between pairs
  of samples, using a choice of distance operators and aggregator
  functions.
}
\usage{
gt.dist(gt1, gt2=gt1, operation='==', aggregator='sum', na.value=NA)
}
\arguments{
  \item{gt1, gt2}{vectors of packed diploid genotypes as from
    \code{\link{fetch.gt.data}}, with 1:1 matching assays.}
  \item{operation}{one of \code{'n'}, \code{'=='}, \code{'!='},
    \code{'ibs==0'}, \code{'ibs==1'}, or \code{'ibs'}.}
  \item{aggregator}{one of \code{'sum'}, \code{'min'}, or \code{'max'}.}
  \item{na.value}{an integer value to use as the result of comparisons
    involving missing genotypes.  The default excludes these
    comparisons.}
}
\details{
  For each genotype assay, an operation is performed on each pairwise
  combination of samples:
  \describe{
    \item{\code{'n'}}{1 if both genotypes are present, else 0.}
    \item{\code{'=='}}{1 if genotypes are identical, else 0.}
    \item{\code{'!='}}{1 if genotypes are not identical, else 0.}
    \item{\code{'ibs==0'}}{1 if genotypes are identical by state (IBS)
      for 0 alleles, else 0.}
    \item{\code{'ibs==1'}}{1 if genotypes are IBS for 1 allele, else 0.}
    \item{\code{'ibs==2'}}{1 if genotypes are IBS for 2 alleles, else 0.}
    \item{\code{'ibs'}}{the number of alleles IBS, 0..2.}
  }
  The matrices of results are then aggregated across genotype assays,
  using one of the available aggregation functions.
}
\value{
  A matrix of integers with \code{length(gt1)} rows and
  \code{length(gt2)} columns, with elements giving the aggregated
  distance for the corresponding pair of samples.
}
\seealso{
  \code{\link{match.gt.data}}, \code{\link{ibd.dataset}}.
}
\examples{
gt <- fetch.gt.data('Demo_2')
gt <- substr(gt$genotype[1:20],1,10)
gt.dist(gt)
gt.dist(gt, operation='ibs', aggregator='min')
}
\keyword{character}
