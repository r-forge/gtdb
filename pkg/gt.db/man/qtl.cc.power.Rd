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
\name{qtl.cc.power}
\alias{qtl.cc.power}
\title{Power Calculation for QTL Case-Control Association Studies}
\description{
  Compute power, or sample size, for a case-control association study,
  with adjustment for linkage disequilibrium, where cases and controls
  are drawn from the tails of a quantitative trait.
}
\usage{
qtl.cc.power(p, add, dom=0, H, upper.tail, lower.tail,
             m, dprime, rsqr, ...)
}
\arguments{
  \item{p}{frequency of the high risk allele.}
  \item{add}{the additive effect per B allele, in units of the
    standard deviation of the within-genotype trait distribution.}
  \item{dom}{the dominance effect: this describes the position of the
    mean trait value for heterozygotes along the interval described by
    the AA and BB homozygotes.  Values range from -1 (purely recessive),
    to 0 (purely additive), to 1 (purely dominant).}
  \item{H}{Fractional heritability attributable to this variant.  This
    is an alternative way of describing the effect size.}
  \item{upper.tail, lower.tail}{the proportions of the population
    distribution of the quantitative trait represented by cases, and
    controls, respectively.}
  \item{m}{if specified, the marker allele frequency.}
  \item{dprime, rsqr}{alternate metrics for specifying linkage
    disequilibrium between the marker and the causal variant.}
  \item{\dots}{additional arguments as in \code{\link{cc.power}}.}
}
\details{
  This is similar in operation to \code{\link{cc.power}}, except that
  cases and controls are defined as tails of the distribution of a
  quantitative trait.
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
\references{
  Schork, N. J., Nath, S. K, Fallin, D., & Chakravarti, A. (2000)
  Linkage disequilibrium analysis of biallelic DNA markers, human
  quantitative trait loci, and threshold-defined case and control
  subjects. \emph{Am. J. Hum. Genet.} \bold{67}: 1208-1218.
}
\seealso{
  \code{\link{qtl.to.cc}}, \code{\link{qtl.power}},
  \code{\link{cc.power}}.
}
\examples{
# reproduce Schork et al. (2000), Table 3, except for factor of 2

mk.table <- function(rows, cols, fn)
{
    .fn <- function(nr,nc)
        do.call(fn, c(cols[nc,,drop=FALSE],rows[nr,,drop=FALSE]))
    x <- lapply(1:nrow(cols), function(nc)
                sapply(1:nrow(rows), .fn, nc))
    x <- unname(do.call('cbind',x))
    colnames(x) <- rownames(cols)
    cbind(rows, x)
}

pwr.fn <- function(p, d, H, dprime, alpha)
    qtl.cc.power(p, dom=d, H=H, upper.tail=0.1, lower.tail=0.1,
                 m=0.25, dprime=dprime, alpha=alpha, power=0.8,
                 method='schork')$N[1]

rows <- expand.grid(dprime=c(0.75,0.5,0.25), H=c(0.1,0.2),
p=c(0.1,0.25))
cols <- expand.grid(alpha=c(0.05,0.0001), d=c(1,-1,0))
mk.table(rows, cols, pwr.fn)
}
\keyword{manip}
