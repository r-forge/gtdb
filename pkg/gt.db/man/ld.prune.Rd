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
\name{ld.prune}
\alias{ld.prune}
\title{Prune SNP List to Limit Linkage Disequilibrium}
\description{
  Takes genotype data as returned by \code{\link{fetch.gt.data}}
  and eliminates SNPs in high LD across a sliding window.
}
\usage{
ld.prune(gt.data, min.maf=0.01, max.rsqr=0.2, span=20, subsets=TRUE)
}
\arguments{
  \item{gt.data}{genotype data for a single contiguous sequence,
    from \code{\link{fetch.gt.data}}.}
  \item{min.maf}{a minimum minor allele frequency for inclusion.}
  \item{max.rsqr}{the maximum allowed R-squared value.}
  \item{span}{the number of sequential SNPs to check for LD.}
  \item{subsets}{an optional list of masks identifying groups of
    subjects within which LD should be checked.}
}
\details{
  LD is computed over a sliding window and SNPs are iteratively
  eliminated to guarantee that each remaining SNP has the specified
  maximum R-squared value for at least \code{span} SNPs on either
  side, ordered by sequence position.

  If multiple \code{subsets} are specified, then the inclusion
  criteria are enforced for ALL subsets.
}
\value{
  A new data frame based on \code{gt.data}, containing a subset of SNPs
  satisfying the LD limitations.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{ld.gt.data}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2')
pt <- fetch.pt.data('Demo_2')
by.plate <- list(pt$plate=='CEU', pt$plate=='YRI',
                 pt$plate=='J+C' & !pt$is_dup)
x <- ld.prune(gt, subsets=by.plate)
nrow(gt)
nrow(x)
}
\keyword{manip}
