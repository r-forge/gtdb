%
% Copyright (C) 2010, 23andMe, Inc.
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
\name{qtl.model}
\alias{qtl.model}
\title{Quantitative Trait Model Parameters}
\description{
  Construct a data structure describing a quantitative trait locus, for
  use in power calculations.
}
\usage{
qtl.model(p, add, dom=0, H)
}
\arguments{
  \item{p}{frequency of the high risk B allele.}
  \item{add}{the additive effect per B allele, in units of the
    standard deviation of the within-genotype trait distribution.}
  \item{dom}{the dominance effect: this describes the position of the
    mean trait value for heterozygotes along the interval described by
    the AA and BB homozygotes.  Values range from -1 (purely recessive),
    to 0 (purely additive), to 1 (purely dominant).}
  \item{H}{Fractional heritability attributable to this variant.  This
    is an alternative way of describing the effect size.}
}
\details{
  This function is used to describe a quantitative trait locus where
  trait values are normally distributed within each genotype, and
  standardized to have a residual variance of 1.  The three genotypes
  have trait distributions centered at \code{-add}, \code{add*dom}, and
  \code{+add}.
}
\value{
  A list with five elements:
  \item{allele.freq}{frequencies of the A and B alleles.}
  \item{gt.freq}{frequencies of AA, AB, BB genotypes.}
  \item{add.effect}{the additive effect per copy of the B allele.}
  \item{dom.effect}{the offset of the mean trait value for the AB
    genotype from the center of the interval between AA and BB
    genotypes.}
  \item{heritability}{the fraction of total variance attributable to
    this QTL.}
}
\seealso{
  \code{\link{qtl.power}}, \code{\link{qtl.to.cc}}.
}
\examples{

}
\keyword{manip}
