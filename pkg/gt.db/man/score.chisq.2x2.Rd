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
\name{score.chisq.2x2}
\alias{score.chisq.2x2}
\title{Chi-Squared Test for Allelic Association}
\description{
  Test for allelic association with a binary outcome, using a
  chi-squared test on a 2x2 table of allele counts.
}
\usage{
score.chisq.2x2(formula, data, ploidy)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}} and details.}
}
\details{
  The 2x2 chi-squared test is the simplest test for association
  between genotypes and a binary outcome.  The model formula should be
  of the form \code{outcome~genotype} without additional terms.  We do
  not apply a continuity correction because that makes the overall
  distribution of P values slightly conservative.

  Sex linked data is handled by counting males as having a single
  allele, which while technically correct may not be reasonable.

  Use of this test is not recommended as it is sensitive to deviations
  from Hardy Weinberg equilibrium.  Usually, \code{\link{score.trend}}
  or \code{\link{score.glm}} are better choices.
}
\value{
  A data frame with one row and three columns:
  \item{pvalue}{P value for the test.}
  \item{effect}{the log odds ratio for the 2x2 table.}
  \item{stderr}{the standard error of the log odds ratio.}
}
\seealso{
  \code{\link{chisq.test}},
  \code{\link{score.chisq}},
  \code{\link{score.trend}},
  \code{\link{score.glm}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~genotype, pt, gt, score.chisq.2x2)
score.gt.data(status~genotype, pt, gt, score.glm)
}
\keyword{htest}

