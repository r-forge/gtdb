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
\name{ibs.gt.data}
\alias{ibs}
\alias{ibs.gt.data}
\alias{ibs.gt.dataset}
\alias{grr}
\title{Calculate Pairwise Identity by State for Genotype Data}
\description{
  Calculates identity-by-state (IBS) for all sample pairs in a set of
  genotype data.
}
\usage{
ibs(x, sample.mask=TRUE)
\method{ibs}{gt.data}(x, sample.mask=TRUE)
\method{ibs}{gt.dataset}(x, sample.mask=TRUE)
grr(x, sample.mask=TRUE)
}
\arguments{
  \item{x}{either a data frame of genotypes from
    \code{\link{fetch.gt.data}}, or a dataset description from
    \code{\link{gt.dataset}}.} 
  \item{sample.mask}{an optional mask identifying a subset of samples
    to be included in the results.}
}
\details{
  \code{ibs.gt.data} and \code{ibs.gt.dataset} assemble counts of
  markers versus the number of alleles IBS.  \code{grr} computes the
  mean and standard deviation of IBS, as described in Abecasis \emph{et
  al.} (2001).
}
\value{
  For \code{ibs.gt.data} and \code{ibs.gt.dataset}, a list of three
  square matrices of counts of markers with IBS=0, IBS=1, and IBS=2, for
  each sample pair.

  For \code{grr}, a list of two square matrices giving the mean and
  standard deviation of the number of alleles IBS, for each sample
  pair.
}
\references{
  Abecasis, G.R., Cherny, Stacey, S.S., Cookson, W.O.C., and Cardon,
  L.R. (2001) GRR: graphical representation of relationship errors.
  \emph{Bioinformatics} \bold{17}: 742-743.
}
\seealso{
  \code{\link{ibd.summary}},
  \code{\link{gt.dataset}}.
}
\examples{
gt.demo.check()
ibs <- ibs(gt.dataset('Demo_1'))
str(ibs)
pt <- fetch.pt.data('Demo_1')
g <- grr(gt.dataset('Demo_1'), (pt$plate=='CEU'))
with(g, xyplot(sd[lower.tri(sd)]~mu[lower.tri(mu)]))
}
\keyword{manip}
