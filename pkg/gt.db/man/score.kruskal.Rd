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
\name{score.kruskal}
\alias{score.kruskal}
\title{Kruskal-Wallis Nonparametric Test for Association}
\description{
  Test for association between genotypes and a quantitative outcome,
  using the nonparametric Kruskal-Wallis test for differences among
  genotype classes.
}
\usage{
score.kruskal(formula, data, ploidy,
              mode=c('general','recessive','dominant'))
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{mode}{genetic mode of action to test.}
}
\details{
  This tests whether the summed ranks of the outcome are independent
  of genotype.  The model formula should be of the form
  \code{outcome~genotype} without additional terms.  The default is to
  make no assumptions about mode of action.
}
\value{
  A data frame with one row and one column:
  \item{pvalue}{asymptotic chi-squared P value for the test.}
}
\seealso{
  \code{\link{kruskal.test}},
  \code{\link{score.jt}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$score <- rnorm(nrow(pt))
score.gt.data(score~genotype, pt, gt, score.kruskal)
score.gt.data(score~genotype, pt, gt, score.lm)
}
\keyword{htest}
