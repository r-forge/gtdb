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
\name{score.coxph}
\alias{score.coxph}
\alias{score.survreg}
\title{Test for Association using Survival Analysis}
\description{
  Test for association using either a semi-parametric Cox proportional
  hazards model, or a parametric survival regression model.
}
\usage{
score.coxph(formula, data, ploidy, test=c('LR','Wald'),
            mode=c('additive','recessive','dominant','general'), ...)
score.survreg(formula, data, ploidy, test=c('LR','Wald'),
              mode=c('additive','recessive','dominant','general'), ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{test}{the test to use for computing P values.}
  \item{mode}{genetic mode of action to test.}
  \item{\dots}{additional arguments passed to
    \code{\link[survival:coxph]{coxph}} or
    \code{\link[survival:survreg]{survreg}}.}
}
\details{
  These functions perform association tests using survival regression
  models.  The model formula may include covariates and should contain a
  single genotype term without interactions.  The left-hand-side of the
  model formula can either identify a survival object, or a numeric
  vector.  If it is numeric, then a survival object will be created
  using the absolute magnitude of the value as the followup time and
  the sign of the value as the event indicator.

  At sex linked loci, haploid males are treated the same as the
  corresponding diploid female homozygotes.
}
\value{
  A data frame with one row and three columns.  An effect size is not
  reported if the mode of action is \code{'general'}.
  \item{pvalue}{P value for the specified \code{test}.}
  \item{effect}{the regression coefficient (log hazard) for the genotype
    term.}
  \item{stderr}{the standard error of the effect size estimate.}
}
\seealso{
  \code{\link[survival:Surv]{Surv}},
  \code{\link[survival:coxph]{coxph}},
  \code{\link[survival:survreg]{survreg}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
library(survival)
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$age <- 50+10*rnorm(nrow(pt))
pt$event <- rbinom(nrow(pt),1,pt$age/100)
pt$survival <- Surv(pt$age,pt$event)
score.gt.data(survival~genotype, pt, gt, score.coxph)
score.gt.data(survival~genotype, pt, gt, score.survreg)
# use numeric coding for outcome
pt$survival <- pt$age * ifelse(pt$event,1,-1)
score.gt.data(survival~genotype, pt, gt, score.survreg)
}
\keyword{htest}
\keyword{regression}
