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
\name{score.glm}
\alias{score.glm}
\title{Test for Association using Logistic Regression}
\description{
  Test for association with a binary outcome using logistic regression.
}
\usage{
score.glm(formula, data, ploidy, drop='genotype',
          test=c('LR','Wald','Rao'),
          mode=c('additive','recessive','dominant','general'))
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{drop}{the model term to test for association.}
  \item{test}{the scoring test.  See details.}
  \item{mode}{genetic mode of action to be tested.}
}
\details{
  The model formula may include covariates and should contain a single
  genotype term without interactions.  For \code{test='LR'} (the default),
  significance is assessed by ANOVA comparing the full model with a null
  model constructed by removing this term.  Other options are
  \code{test='Wald'} (a Wald test on the coefficient of the genotype
  term in the regression model), and \code{test='Rao'} (the Rao score
  test, which is equivalent to the Cochran-Armitage trend test if there
  are no covariates).

  At sex linked loci, haploid males are treated the same as the
  corresponding diploid female homozygotes.
}
\value{
  A data frame with one row and three columns.  An effect size is not
  reported if the mode of action is \code{'general'}.
  \item{pvalue}{P value for an F test comparing the full and null models.}
  \item{effect}{the regression coefficient (log odds) for \code{term}.}
  \item{stderr}{the standard error of the effect size estimate.}
}
\seealso{
  \code{\link{lm}},
  \code{\link{score.trend}},
  \code{\link{score.chisq}},
  \code{\link{score.glm.general}},
  \code{\link{score.glm.groups}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~plate+genotype, pt, gt, score.glm)
score.gt.data(status~plate+genotype, pt, gt, score.glm.general)
}
\keyword{htest}
\keyword{regression}
