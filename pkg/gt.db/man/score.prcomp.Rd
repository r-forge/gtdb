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
\name{score.prcomp}
\alias{score.prcomp}
\title{Test Phenotypic Association with Principal Components}
\description{
  Evaluates evidence for association between sample loadings from a
  principal components analysis, and one or more phenotypes.
}
\usage{
score.prcomp(formula, pc.data, pt.data, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fitted.
    The left hand side should be given as \code{PC}.}
  \item{pc.data}{a PCA result structure returned by
    \code{\link{prcomp.gt.data}} or \code{\link{prcomp.gt.dataset}}.}
  \item{pt.data}{phenotype data from \code{\link{fetch.pt.data}}.}
  \item{\dots}{additional arguments passed to \code{\link{lm}}.}
}
\details{
  The specified model formula is evaluated by linear regression for
  each principal component, and results are summarized as in an ANOVA
  based on the proportion of variance of the principal component
  explained.
}
\value{
  A data frame with one row per principal component and three columns:
  \item{R\^2}{the adjusted multivariate R\^2 statistic for the model.}
  \item{F value}{an F statistic for the model.}
  \item{Pr(>F)}{the tail probability associated with this F statistic.}
}
\seealso{
  \code{\link{prcomp.gt.data}}, \code{\link{prcomp.gt.dataset}},
  \code{\link{lm}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_1')
pt <- fetch.pt.data('Demo_1')
pc <- prcomp(subset(gt, (ploidy=='A')),
             sample.mask=(is.na(pt$father)))
score.prcomp(PC~panel, pc, pt)
}
\keyword{manip}
