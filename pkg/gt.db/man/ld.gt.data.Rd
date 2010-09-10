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
\name{ld.gt.data}
\alias{ld.gt.data}
\title{Compute Pairwise Linkage Disequilibrium}
\description{
  Takes genotype data as returned by \code{\link{fetch.gt.data}}
  and computes a matrix of pairwise linkage disequilibrium values.
}
\usage{
ld.gt.data(g1, g2=g1, outer=TRUE,
           measure=c('rsqr','dprime','delta','pvalue',
                     'lod','none','failed'),
           method=c('iterate','exact'), epsilon=1e-6, max.it=100)
}
\arguments{
  \item{g1, g2}{genotype data from \code{\link{fetch.gt.data}}.}
  \item{outer}{logical: specifies whether result should be an outer
    product of \code{g1} and \code{g2}.}
  \item{measure}{the measure of LD to report.  See details.}
  \item{method}{the method to use for calculating LD.}
  \item{max.it}{the maximum number of EM iterations to perform.}
  \item{epsilon}{convergence criterion for the EM algorithm.}
}
\details{
  LD is computed between unphased genotypes in \code{g1} and \code{g2}
  under an assumption of Hardy Weinberg equilibrium.  Two algorithms are
  implemented: an iterative EM solution, and an exact solution.

  All markers should have consistent ploidy status.  The LD
  calculations account for haploid and diploid genotypes by using
  gender information at sex-linked loci.

  With the iterative algorithm, a warning is printed if some pairwise LD
  calculations do not converge in the specified number of iterations.
  If \code{measure='failed'}, then the function returns an array of
  logical values indicating which cases failed to converge.

  For the exact algorithm, the most-likely result is reported.  With
  extreme deviations from HWE, it is possible for more than one
  solution to have the same likelihood, and in these cases, the one
  that is reported is essentially arbitrary.
}
\value{
  If \code{outer} is \code{TRUE}, then the result is a matrix of LD
  values with one row per element of \code{g1} and one column per
  element of \code{g2}.  If \code{outer} is \code{FALSE}, then the
  result is a vector of values obtained by comparing corresponding
  elements of \code{g1} with \code{g2} with recycling of the shorter
  argument.  This will usually only make sense when the lengths of
  \code{g1} and \code{g2} are either equal or equal to \code{1}.

  The number of pairwise calculations that failed to converge is
  returned in an attribute with the name \code{'failed'}.
}
\references{
  Weir, B.S. (1996) \emph{Genetic Data Analysis II}.  Sinauer,
  Sunderland, MA.

  Gaunt, T. R., Rodriguez, S., & Day, I. N. M. (2007) Cubic exact
  solutions for the estimation of pairwise haplotype frequencies:
  implications for linkage disequilibrium analyses and a web tool
  'CubeX'.  \emph{BMC Bioinformatics} \bold{8}: 428.

  CubeX. \url{http://www.oege.org/software/cubex/}.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{ld.plot}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2')[1:6,]
ld.gt.data(gt, measure='rsqr')
ld.gt.data(gt, measure='dprime')
ld.gt.data(gt, measure='pvalue')
ld.gt.data(gt[1:3,], gt[3:5,], outer=TRUE)
ld.gt.data(gt[1:3,], gt[3:5,], outer=FALSE)
}
\keyword{manip}
