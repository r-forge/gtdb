%
% Copyright (C) 2010, 23andMe, Inc.
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
\name{qtl.to.cc}
\alias{qtl.to.cc}
\title{Create Case-Control Model from Quantitative Trait Tails}
\description{
  Construct a binary outcome model where cases and controls are drawn
  from the upper and lower tails of the distribution of a quantitative
  trait.
}
\usage{
qtl.to.cc(model, upper.tail, lower.tail, pdist=pnorm)
}
\arguments{
  \item{model}{a model from \code{\link{qtl.model}}.}
  \item{upper.tail, lower.tail}{the proportions of the population
    distribution of the quantitative trait represented by cases, and
    controls, respectively.}
  \item{pdist}{the cumulative probability distribution of the trait
    within each genotype.}
}
\details{
}
\value{
  A list consisting of the original quantitative trait model, with
  additional elements as returned by \code{\link{cc.model}} describing
  the derived binary outcome.
}
\references{
  Schork, N. J., Nath, S. K, Fallin, D., & Chakravarti, A. (2000)
  Linkage disequilibrium analysis of biallelic DNA markers, human
  quantitative trait loci, and threshold-defined case and control
  subjects. \emph{Am. J. Hum. Genet.} \bold{67}: 1208-1218.
}
\seealso{
  \code{\link{qtl.cc.power}}, \code{\link{qtl.model}},
  \code{\link{cc.model}}.
}
\examples{

}
\keyword{manip}
