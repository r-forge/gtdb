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
\name{qtl.power}
\alias{qtl.power}
\title{Power Calculation for Quantitative Trait Association Studies}
\description{
  Compute power, or sample size, for detecting association with a
  quantitative trait, where genotype only affects the location of the
  trait distribution.
}
\usage{
qtl.power(p, add, dom=0, H, N, alpha=0.05,
          method=c('model','simulate'), ...)
}
\arguments{
  \item{p}{frequency of the high risk allele.}
  \item{add}{the additive effect per allele, in units of the
    standard deviation of the within-genotype trait distribution.}
  \item{dom}{the dominance effect: this describes the position of the
    mean trait value for heterozygotes along the interval described by
    the two homozygotes.  Values range from -1 (purely recessive), to
    0 (purely additive), to 1 (purely dominant).}
  \item{H}{Fractional heritability attributable to this variant.  This
    is an alternative way of describing the effect size.}
  \item{N}{The total number of individuals to be genotyped.}
  \item{alpha}{the desired significance level.}
  \item{power}{if specified, the desired study power.}
  \item{method}{the scoring criterion used to determine power.  See
    details.}
  \item{\dots}{additional method-specific arguments.  See details.}
}
\details{
  This calculates study power for a given sample size (if \code{N} is
  specified), or sample size for a desired power (if \code{power} is
  specified).

  Two methods for calculating power are implemented:
  \describe{
    \item{model}{model-based power calculation for detecting an effect
      by linear regression.}
    \item{simulate}{simulated power for an arbitrary scoring function.}
  }

  The \code{simulate} method has the following additional arguments:
  \describe{
    \item{score.fn}{the scoring function to use: defaults to
      \code{\link{score.lm}}.}
    \item{tries}{the number of iterations to perform.}
    \item{progress}{logical: indicates whether to show a progress bar.}
    \item{rdist}{function for generating the residual within-genotype
      error, defaults to \code{\link{rnorm}}.  If specified, this should
      have a variance of 1.}
    \item{\dots}{additional arguments passed to \code{score.fn}}
  }

}
\value{
  A list with four elements:
  \item{model}{a list as returned by \code{\link{cc.model}}.}
  \item{N}{the study size (numbers of cases and controls).}
  \item{alpha}{the desired false positive rate (probability of
    incorrectly rejecting the null hypothesis when it is true).}
  \item{power}{the study power (probability of rejecting the null
    hypothesis).}
}
\seealso{
  \code{\link{qtl.model}}, \code{\link{score.lm}},
  \code{\link{qtl.cc.power}}.
}
\examples{
# power to detect association comparing upper and lower tails
sapply(seq(0.05,0.5,0.05), function(x)
       qtl.cc.power(p=0.1, add=0.2, N=2000*x,
                    upper.tail=x, lower.tail=x)$power)
# power for detecting QTL by linear regression
qtl.power(p=0.1, add=0.2, N=1000)$power
# power for non-parametric tests
qtl.power(p=0.1, add=0.2, N=1000, method='simulate',
          score.fn=score.kruskal, tries=200)$power
qtl.power(p=0.1, add=0.2, N=1000, method='simulate',
          score.fn=score.jt, tries=200)$power
}
\keyword{manip}
