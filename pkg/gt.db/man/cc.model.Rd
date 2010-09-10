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
\name{cc.model}
\alias{cc.model}
\title{Case-Control Disease Model Parameters}
\description{
  Construct a data structure describing a dichotomous disease
  phenotype, for use in power calculations.
}
\usage{
cc.model(p, prevalence, low.risk, rel.risk,
         odds.ratio, pop.controls=FALSE)
}
\arguments{
  \item{p}{frequency of the high risk B allele.}
  \item{prevalence}{population prevalence of case status.}
  \item{low.risk}{baseline risk for the AA genotype.}
  \item{rel.risk}{either the per-allele (log-additive) relative risk,
    or a vector of two relative risks for AB vs AA, and BB vs AA.}
  \item{odds.ratio}{either the per-allele (log-additive) odds ratio,
    or a vector of two odds ratios for AB vs AA, and BB vs AA.}
  \item{pop.controls}{logical: indicates if controls are unverified and
    are assumed to include cases at the population prevalence.}
}
\details{
  To describe the distribution of risk versus genotype, the caller
  needs to specify either \code{prevalence} or \code{low.risk}, and
  either \code{rel.risk} or \code{odds.ratio}.
}
\value{
  A list with the following elements:
  \item{prevalence}{population prevalence of case, control status.}
  \item{allele.freq}{frequencies of the A and B alleles.}
  \item{gt.freq}{frequencies of AA, AB, BB genotypes.}
  \item{penetrance}{penetrance for AA, AB, BB genotypes for case,
    control status.}
  \item{odds.ratio}{odds ratios for AB vs AA, BB vs AA.}
  \item{rel.risk}{genotype relative risks for AB vs AA, BB vs AA.}
  \item{exp.freq}{expected allele frequencies in cases, controls.}
  \item{exp.gt.freq}{expected genotype frequencies in cases, controls.}
}
\seealso{
  \code{\link{cc.power}}, \code{\link{adjust.ld}}.
}
\examples{
cc.model(0.1, prevalence=0.2, odds.ratio=1.5)
}
\keyword{manip}
