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
\name{score.trend}
\alias{score.trend}
\title{Cochran-Armitage Trend Test for Association}
\description{
  Test for genotype assocation with a binary outcome, using the
  Cochran-Armitage test for trend in proportions.
}
\usage{
score.trend(formula, data, ploidy)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
}
\details{
  The Cochran-Armitage test for trend is a test for a linear
  dependence of outcome proportions on genotype expressed as an allele
  count.  The model formula should be of the form
  \code{outcome~genotype} without additional terms.

  In place of a conventional effect size, we return the slope and
  standard error for the trend in proportions.  This slope is thus the
  fitted change in outcome proportion per allele.
}
\value{
  A data frame with one row and three columns:
  \item{pvalue}{P value for test.}
  \item{effect}{the slope of the trend in proportions.}
  \item{stderr}{the standard error of the slope.}
}
\seealso{
  \code{\link{prop.trend.test}},
  \code{\link{score.glm}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~genotype, pt, gt, score.trend)
score.gt.data(status~genotype, pt, gt, score.glm)
}
\keyword{htest}

