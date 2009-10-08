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
\name{ld.plot}
\alias{ld.plot}
\title{Pairwise Linkage Disequilibrium Plot}
\description{
  Generates a level plot of the strength of pairwise linkage
  disequilibrium among a set of SNPs.
}
\usage{
ld.plot(gt.data, col=gray(seq(1,0,-0.01)), measure='rsqr', 
        rotate=FALSE, equal=TRUE, colorkey=list(height=0.5),
        scales, ...)
}
\arguments{
  \item{gt.data}{genotype data from \code{\link{fetch.gt.data}}.}
  \item{col}{a color gradient.}
  \item{measure}{the measure(s) of LD to be plotted.  See details.}
  \item{rotate}{logical: indicates if the plot should be rotated
    clockwise 45 degrees.}
  \item{equal}{logical: indicates if the LD matrix should be displayed
    with equal SNP spacing, or in genomic coordinates.}
  \item{colorkey, scales}{passed to
    \code{\link[lattice:levelplot]{levelplot}}.}
  \item{\dots}{additional arguments to pass to
    \code{\link[lattice:levelplot]{levelplot}}.}
}
\details{
SNPs are sorted by accession and contig position.  Values of
\code{measure} are as in \code{\link{ld.gt.data}}.  If \code{measure}
has two elements, then the first value is plotted in the upper left
triangle and the second value is plotted in the lower right triangle.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{ld.gt.data}} \code{\link[lattice:levelplot]{levelplot}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2')
pt <- fetch.pt.data('Demo_2')
gt <- mask.gt.data(gt, pt$plate=='CEU')
gs <- summary.gt.data(gt)
ld.plot(subset(gt, pass & gs$freq.a>0.05 & gs$freq.b>0.05),
        measure=c('dprime','rsqr'), rotate=TRUE, equal=FALSE)
}
\keyword{hplot}
