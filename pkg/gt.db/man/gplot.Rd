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
\name{gplot}
\alias{gplot}
\title{Genome-Wide Level Plot}
\description{
  Generates a level plot of genome-wide data organized by chromosome
  and position.
}
\usage{
gplot(formula, data, aggr.fn=max, rescale=FALSE, binsz=1e6,
      subset=TRUE, col.regions=rev(heat.colors(100)[10:100]),
      scales=list(x=list(at=seq(0,250,20), draw=TRUE)),
      shrink=list(x=1,y=0.75), colorkey=list(height=0.25),
      xlab='Position, Mb', ylab='Chromosome', ...)
}
\arguments{
  \item{formula}{a one-sided formula (i.e. \code{~ x}) describing the
                 values to be plotted.}
  \item{data}{a data frame containing values required to evaluate
              \code{formula}, as well as \code{chromosome} and 
              \code{position} columns.}
  \item{aggr.fn}{an aggregating function to apply to data within each bin.}
  \item{rescale}{logical: indicates if results should be scaled based
                 on the genome-wide average bin value.}
  \item{binsz}{bin size in base pairs}
  \item{subset}{a logical or integer vector identifying rows of 
                \code{data} to be included in the plot.}
  \item{col.regions, scales, shrink, colorkey, xlab, ylab}{see
    \code{\link[lattice:levelplot]{levelplot}}.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:levelplot]{levelplot}}.}
}
\details{
  The supplied formula is first evaluated across the input data, and 
  then the aggregating function \code{aggr.fn} is applied to results
  in bins of genomic coordinates.  Appropriate aggregating functions
  include \code{\link{max}}, \code{\link{min}}, \code{\link{sum}}, etc.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link[lattice:levelplot]{levelplot}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_1')
gt <- fetch.gt.data('Demo_1')
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
r <- score.gt.data(status~genotype, pt, gt, score.chisq)
gplot(~-log10(pvalue), merge(r,gt))
}
\keyword{hplot}
