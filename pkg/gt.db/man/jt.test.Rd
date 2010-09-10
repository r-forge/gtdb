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
\name{jt.test}
\alias{jt.test}
\title{Jonckheere-Terpstra Nonparametric Test for Trend}
\description{
  Test for association between genotypes and a quantitative outcome,
  using the nonparametric Jonckheere-Terpstra test for ordered
  differences among genotype classes.
}
\usage{
jt.test(x, y, alternative=c("two.sided", "decreasing", "increasing"),
        asymp=TRUE, correct=FALSE, perm=0, na.action=c("omit","fail"),
        permgraph=FALSE, permreps=FALSE)
}
\arguments{
  \item{x}{vector of quantitative response values.}
  \item{y}{group membership.}
  \item{alternative}{alternative hypothesis to be tested.}
  \item{asymp}{logical: use asymptotic formula for variance, or don't bother.}
  \item{correct}{logical: apply continuity correction?}
  \item{perm}{number of repetitions for a permutation test.}
  \item{na.action}{what to do with missing data.}
  \item{permgraph}{logical: draw a histogram of the permutations?}
  \item{permreps}{logical: return the permutations?}
}
\details{
  The Jonckheere-Terpstra test computes Mann-Whitney rank sum
  statistics for ordered group labels, and does a two-sided test on
  the sum of those statistics.  This tests for a monotonic trend in
  response as a function of group membership.  The model formula
  should be of the form \code{outcome~genotype} without additional
  terms.
}
\value{
  A list of class htest, with the following components:
  \item{statistic}{the observed J-T statistic.}
  \item{alternative}{same as input.}
  \item{method}{the string: \code{"Jonckheere-Terpstra test"}.}
  \item{data.name}{the names of the input data.}
  \item{EH}{the expected test statistic based on sample size.}
  \item{VH}{variance (adjusted for ties if necessary).}
  \item{p.value}{asymptotic p-value.}
}
\references{
  \url{http://tolstoy.newcastle.edu.au/R/help/06/06/30112.html}.
}
\seealso{
  \code{\link{kruskal.test}},
  \code{\link{score.kruskal}}.
}
\examples{
x <- rnorm(30) + c(rep(1,10),rep(2,20))
y <- c(rep(1,10),rep(2,10),rep(3,10))
jt.test(x, y)
}
\keyword{htest}

