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
\name{manhattan.plot}
\alias{manhattan.plot}
\title{Genome-Wide Manhattan Plot}
\description{
  Generates a "Manhattan plot" of genome-wide data, with chromosomes
  arranged in ascending order along the X axis.
}
\usage{
manhattan.plot(y, data, gap=0, threshold=-log10(5e-8), around=0,
               xticks=c(1:12,14,16,18,20,22,'X','Y'), cex=0.25,
               xlab=, ylab=expression(-log[10](pvalue)), yrange, ylim,
               scale=c('Mb','Kb','bp'),
               col=c('#d0d0d0','#e0e0e0','#ff0000'), ...)
}
\arguments{
  \item{y}{the scores, evaluated in context of \code{data}.}
  \item{data}{a data frame that includes \code{scaffold} and
    \code{position} columns giving genomic positions.}
  \item{gap}{a gap size, in base pairs, to insert between chromosomes.}
  \item{threshold}{a score threshold for a \dQuote{hit}.}
  \item{around}{in base pairs, specifies that SNPs within this interval
    of a hit should also be highlighted.}
  \item{xticks}{a list of chromosomes to be labeled in the plot.}
  \item{cex, xlab, ylab, ylim}{see \code{\link[lattice:xyplot]{xyplot}}.}
  \item{yrange}{similar to \code{ylim}, except that this range will be
    padded as if it were the range of actual data points.}
  \item{scale}{determines the scaling of tick labels on the X axis, if
    the data all falls on a single chromosome.}
  \item{col}{colors to use for odd and even chromosomes, and hits.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:xyplot]{xyplot}}.}
}
\details{
  Typically, the scores are expected to be \code{-log10(pvalue)}.  If
  the source data spans multiple chromosomes, then chromosome names are
  shown along the X axis.  If the data is all on a single chromosome,
  then chromosomal positions are shown.

  The default scaling of the Y axis if neither \code{yrange} or
  \code{ylim} are specified will place the threshold level at least
  10% below the upper edge of the plot area.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{gplot}}, \code{\link[lattice:xyplot]{xyplot}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_1')
gt <- fetch.gt.data('Demo_1')
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
r <- score.gt.data(status~genotype, pt, gt, score.chisq)
manhattan.plot(-log10(pvalue), merge(r,gt), threshold=-log10(0.01))
}
\keyword{hplot}
