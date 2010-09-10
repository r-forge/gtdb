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
\name{score.jt}
\alias{score.jt}
\title{Jonckheere-Terpstra Nonparametric Test for Association}
\description{
  Test for association between genotypes and a quantitative outcome,
  using the nonparametric Jonckheere-Terpstra test for ordered
  differences among genotype classes.
}
\usage{
score.jt(formula, data, ploidy, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{\dots}{additional arguments passed to \code{\link{jt.test}}.}
}
\details{
  The Jonckheere-Terpstra test computes Mann-Whitney rank sum
  statistics for ordered diploid genotypes, and does a two-sided test
  on the sum of those statistics.  This tests for a monotonic trend in
  outcomes as a function of genotype.  The model formula should be of
  the form \code{outcome~genotype} without additional terms.

  The test does not directly provide an estimate of an effect size.
  Instead, we center and rescale the test statistic to fall in the
  range \code{-1} to \code{+1}.  The expected value of this score is
  independent of sample size and indicates the direction of effect.
}
\value{
  A data frame with one row and three columns:
  \item{pvalue}{P value for the test assuming an asymptotic normal
    distribution for the test statistic.}
  \item{effect}{the test statistic centered and scaled to be in the
    range \code{-1} to \code{+1}.}
  \item{stderr}{the estimated standard error of \code{effect}.}
}
\references{
  \url{http://tolstoy.newcastle.edu.au/R/help/06/06/30112.html}.
}
\seealso{
  \code{\link{jt.test}},
  \code{\link{score.kruskal}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$score <- rnorm(nrow(pt))
score.gt.data(score~genotype, pt, gt, score.kruskal)
score.gt.data(score~genotype, pt, gt, score.jt)
}
\keyword{htest}

