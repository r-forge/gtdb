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
\name{score.glm.groups}
\alias{score.glm.groups}
\title{Test for Association using a Logistic Model with Subgroup Effects}
\description{
  Test for association with a binary outcome, using logistic
  regression assuming a log-additive allelic effect on risk for the
  genotype term, where effects are allowed to vary by subgroup.
}
\usage{
score.glm.groups(formula, data, ploidy, quick=FALSE,
                 mode=c('additive','recessive','dominant'))
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{ploidy}{see \code{\link{fetch.gt.data}}.  Ignored.}
  \item{quick}{logical: if set, just return the overall test result.}
  \item{mode}{genetic mode of action to test.}
}
\details{
  The model formula can contain multiple genotype terms with
  interactions, where the interacting variables are factors.  Overall
  significance is assessed by a likelihood ratio test comparing the
  full model with a null model constructed by removing all genotype
  terms.  Significance of the interaction term(s) is assessed by a
  likelihood ratio test comparing the full model with a model with no
  interactions.  Finally, genotype effects are computed for each
  subgroup.  At sex linked loci, haploid males are treated the same as
  the corresponding diploid female homozygotes.
}
\value{
  A data frame with four columns and at least two rows.  The first two
  rows describe likelihood ratio tests for the full model versus a
  null model, and the full model versus a model with no subgroup
  interaction.  Additional rows describe tests for effects within each
  subgroup.
  \item{term}{for subgroup effects, the corresponding
    genotype:subgroup interaction term.}
  \item{pvalue}{for tests on nested models, the P value for a
    likelihood ratio test comparing the two models.  For subgroup
    effects, the P value of a Wald test on the effects size.}
  \item{effect}{for subgroups, the estimated allelic effect (log odds)
    from the regression.}
  \item{stderr}{for subgroups, the estimated standard error of the
    effect size.} 
}
\seealso{
  \code{\link{glm}},
  \code{\link{score.glm}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:5,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~plate+genotype, pt, gt, score.glm)
score.gt.data(status~plate*genotype, pt, gt, score.glm.groups)
}
\keyword{htest}
\keyword{regression}
