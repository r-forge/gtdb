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
\name{score.lm}
\alias{score.lm}
\title{Test for Association using a Simple Linear Model}
\description{
  Test for association using linear regression.
}
\usage{
score.lm(formula, data, ploidy, drop='genotype',
         mode=c('additive','recessive','dominant','general'))
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{drop}{the model term to test for association.}
  \item{mode}{genetic mode of action to test.}
}
\details{
  The model formula may include covariates and should contain a single
  genotype term without interactions.  Significance is assessed by
  ANOVA comparing the full model with a null model constructed by
  removing this term. 

  At sex linked loci, haploid males are treated the same as the
  corresponding diploid female homozygotes.
}
\value{
  A data frame with one row and three columns.  An effect size is not
  reported if the mode of action is \code{'general'}.
  \item{pvalue}{P value for an F test comparing the full and null models.}
  \item{effect}{the regression coefficient for \code{term}.}
  \item{stderr}{the standard error of the effect size estimate.}
}
\seealso{
  \code{\link{lm}},
  \code{\link{score.lm.general}},
  \code{\link{score.lm.groups}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$score <- rnorm(nrow(pt))
score.gt.data(score~plate+genotype, pt, gt, score.lm)
score.gt.data(score~genotype, pt, gt, score.lm, mode='recessive')
score.gt.data(score~genotype, pt, gt, score.lm, mode='general')
score.gt.data(score~genotype, pt, gt, score.lm.general)
}
\keyword{htest}
\keyword{regression}
