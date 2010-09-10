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
\name{score.gt.data}
\alias{score.gt.data}
\title{Test SNPs for Association}
\description{
  Perform a series of single-point SNP association tests, using an
  arbitrary scoring function.
}
\usage{
score.gt.data(formula, pt.data, gt.data, score.fn=NULL,
              pt.filter=TRUE, gt.filter=TRUE, dosage=FALSE,
              progress=FALSE, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{pt.data}{a data frame of phenotypes from \code{\link{fetch.pt.data}}.}
  \item{gt.data}{a data frame of genotypes from \code{\link{fetch.gt.data}}.}
  \item{score.fn}{the scoring function to be applied to each SNP.}
  \item{pt.filter}{an expression to use for subsetting on the
    phenotype table.}
  \item{gt.filter}{an expression to use for subsetting on the
    genotype table.}
  \item{dosage}{logical: specifies whether to use allele dosage
    in association tests in place of genotype scores.}
  \item{progress}{logical: specifies if a progress bar should be
    displayed on the R console.}
  \item{\dots}{additional arguments to pass to \code{score.fn}.}
}
\details{
  For each row of genotypes in \code{gt.data}, a data frame will be
  constructed by merging columns of \code{pt.data} referenced in the
  model formula with those genotypes.  The specified score function
  will be invoked for each of these data frames in turn.

  If no scoring function is specified, then one will be deduced from
  the form of the model formula.  A trend test will be used for simple
  binary-outcome models of the form "status~genotype".  Logistic
  regression will be used for binary outcomes with more complicated
  model functions.  And linear regression will be used for
  quantitative outcomes.

  The filters \code{pt.filter} and \code{gt.filter} are evaluated
  in the context of the \code{pt.data} and \code{gt.data} tables,
  respectively, similar to how \code{subset} works.

  It is expected that some tests will fail for some SNPs (say, due to
  degeneracies where a SNP or trait is invariant across the samples
  for which data is available).  Errors from the scoring function are
  converted to warnings.
}
\value{
  A data frame with up to five columns, and one or more rows per SNP
  tested, depending on the form of the test.  Some columns may not be
  populated for certain test types.
  \item{assay.name}{an identifier for the assay tested.}
  \item{term}{for tests that return more than one result per assay,
    something specifying the individual results.}
  \item{pvalue}{an uncorrected P value for this test.}
  \item{effect}{an estimated effect size, where that makes sense.}
  \item{stderr}{the estimated standard error of the effect size.}
}
\seealso{
  \code{\link{score.chisq.2x2}},
  \code{\link{score.chisq}},
  \code{\link{score.fisher}},
  \code{\link{score.trend}},
  \code{\link{score.kruskal}},
  \code{\link{score.jt}},
  \code{\link{score.lm}},
  \code{\link{score.lm.general}},
  \code{\link{score.lm.groups}},
  \code{\link{score.glm}},
  \code{\link{score.glm.general}},
  \code{\link{score.glm.groups}},
  \code{\link{score.and.store}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')[1:10,]
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
score.gt.data(status~genotype, pt, gt)
score.gt.data(status~plate+genotype, pt, gt, score.glm,
              pt.filter=(is.na(father)))
}
\keyword{htest}
\keyword{regression}
