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
\name{score.glm.general}
\alias{score.glm.general}
\title{Association Test with a Logistic Model and General Mode of Action}
\description{
  Test for association using logistic regression, with a general mode
  of action.
}
\usage{
score.glm.general(formula, data, ploidy)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
}
\details{
  This is similar to \code{\link{score.glm}}, except that it considers
  a general mode of action with (log) additive and dominance effects.
  Overall significance is still assessed by ANOVA comparing the full
  model with a null model constructed by removing this term. 
}
\value{
  A data frame with four columns.  If there are only two distinct
  genotypes, then results are equivalent to \code{score.glm} with
  additive mode of action.  If there are three distinct genotypes,
  then the first row of results describes the overall ANOVA test;
  the second describes the fitted additive effect; and the third
  describes the fitted dominance effect.
  \item{term}{\code{NA}, \code{'additive'}, or \code{'dominance'}.}
  \item{pvalue}{For the overall test, the ANOVA test result.
    For additive and dominance effects, Wald test results.}
  \item{effect}{the regression coefficient, for the additive or
    dominant effect.}
  \item{stderr}{the standard error of the effect size estimate.}
}
\seealso{
  \code{\link{glm}},
  \code{\link{score.glm}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~genotype, pt, gt, score.glm)
score.gt.data(status~genotype, pt, gt, score.glm.general)
}
\keyword{htest}
\keyword{regression}
