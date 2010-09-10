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
\name{score.fisher}
\alias{score.fisher}
\title{Fisher's Exact Test for Genotypic Association}
\description{
  Performs Fisher's exact test for genotype association with a
  categorical outcome, on an NxM table of genotype counts.
}
\usage{
score.fisher(formula, data, ploidy,
             mode=c('general','recessive','dominant'), ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{mode}{genetic mode of action to test.}
  \item{\dots}{additional arguments passed to \code{\link{fisher.test}}.}
}
\details{
  This evaluates Fisher's exact test for association on a contingency
  table for a categorical outcome (with any number of levels) versus
  diploid genotypes.  The model formula should be of the form
  \code{outcome~genotype} without additional terms.  Empty levels for
  the outcome or the genotype will be dropped.

  Sex linked data is handled by treating males as diploid homozygotes
  for the corresponding haploid alleles, which may or may not be a
  reasonable thing to do.
}
\value{
  A data frame with one row and two columns:
  \item{pvalue}{P value for the test.}
  \item{effect}{for 2x2 tables, the log odds ratio.}
}
\seealso{
  \code{\link{fisher.test}},
  \code{\link{score.chisq}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~genotype, pt, gt, score.chisq)
score.gt.data(status~genotype, pt, gt, score.fisher)
}
\keyword{htest}
