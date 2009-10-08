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
\name{hwe.test}
\alias{hwe.test}
\title{Tests for Hardy Weinberg Equilibrium}
\description{
  Asymtotic and exact tests for Hardy Weinberg equilibrium conditional
  on observed marginal allele frequencies, for biallelic genotype data.
}
\usage{
hwe.test(aa, ab, bb, test=c('lratio','chisq','exact'),
         tail=c('both','lower','upper'))
}
\arguments{
  \item{aa, ab, bb}{biallelic genotype counts (may be vectors).}
  \item{test}{the type of test to perform.  See details.}
  \item{tail}{See details.}
}
\details{
  Three tests are implemented: the likelihood ratio test, the
  traditional chi-squared test without continuity correction, and the
  conditional exact test.  All are defined as in Weir (1996).  All the
  tests are vectorized for efficiency.

  The default is to perform a two-sided test.  If \code{tail='lower'},
  then the result is the probability of seeing no more than the
  observed number of heterozygotes under the hypothesis of Hardy
  Weinberg equilibrium.  If \code{tail='upper'}, then the result is
  the probability of seeing at least as many as the observed number of
  heterozygotes.
}
\value{
  A vector of P values.
}
\references{
  Weir, B.S. (1996) \emph{Genetic Data Analysis II}.
  Sinauer, Sunderland, MA.
}
\examples{
# Table 3.1 in Weir (1996)
d <- data.frame(aa=0:9, ab=seq(19,1,-2), bb=21:30)
d$p.exact <- with(d, hwe.test(aa,ab,bb,'exact'))
d$p.chisq <- with(d, hwe.test(aa,ab,bb,'chisq'))
d$chisq <- qchisq(d$p.chisq,df=1,lower.tail=FALSE)
d[order(d$p.chisq),]
#
gt.demo.check()
gt <- fetch.gt.data('Demo_1')
pt <- fetch.pt.data('Demo_1')
s <- summary.gt.data(gt, pt$plate=='J+C')
d <- with(s,data.frame(
    lratio=hwe.test(AA,AB,BB,'lratio'),
    chisq=hwe.test(AA,AB,BB,'chisq'),
    exact=hwe.test(AA,AB,BB,'exact')
))
p <- list(limits=c(-10,0), at=c(-2,-4,-6,-8),
          labels=c('1E-2','1E-4','1E-6','1E-8'))
splom(~log(d,10), pscales=list(p,p,p))
}
\keyword{htest}
