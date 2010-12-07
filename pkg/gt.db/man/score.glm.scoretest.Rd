%
% Copyright (C) 2010, 23andMe, Inc.
% 
% Written by David A. Hinds <dhinds@23andMe.com>
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
\name{score.glm.scoretest}
\alias{score.glm.scoretest}
\title{Fast Score Test for Association using Logistic Regression}
\description{
  A score test for association with a binary outcome based on logistic
  regression, using a precalculated null model.
}
\usage{
score.glm.scoretest(formula, data, ploidy,
                    precalc=glm(.null.model(formula),binomial,data))
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{precalc}{a null model to be precalculated once.}
}
\details{
  This is a variation on \code{\link{score.glm}} that calculates P
  values based on Rao's score test.  It precalculates the base model
  (which can contain covariates).  In order to use this precalculated
  model, there can be no missing genotype data, so missing observations
  are replaced by the mean.  A consequence is that only the additive
  mode of action is supported.  It is much faster than
  \code{\link{score.glm}}.
  
  At sex linked loci, haploid males are treated the same as the
  corresponding diploid female homozygotes.
}
\value{
  A data frame with one row and three columns:
  \item{pvalue}{P value for the score test.}
  \item{effect}{always \code{NA}.}
  \item{stderr}{always \code{NA}.}
}
\seealso{
  \code{\link{glm}}, \code{\link[statmod:glm.scoretest]{glm.scoretest}},
  \code{\link{score.glm}}, \code{\link{score.trend}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~genotype, pt, gt, score.glm, test='Rao')
score.gt.data(status~genotype, pt, gt, score.trend)
score.gt.data(status~genotype, pt, gt, score.glm.scoretest)
}
\keyword{htest}
\keyword{regression}
