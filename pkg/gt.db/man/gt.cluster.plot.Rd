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
\name{gt.cluster.plot}
\alias{gt.cluster.plot}
\alias{xyplot.gt.data}
\title{Plot Genotype Cluster Data}
\description{
  Plot raw genotype cluster data (i.e. signal intensities), with
  predicted ellipsoid boundary regions.
}
\usage{
gt.cluster.plot(data, by=assay.name, rescale=TRUE, bounds=c(0.5,0.95),
                min.points=4, between=list(x=0.5,y=0.5),
                scales=list(alternating=0), xlab=NULL, ylab=NULL,
                par.settings=.gt.settings, ...)
\method{xyplot}{gt.data}(x, ...)
}
\arguments{
  \item{data}{an unpacked data frame of genotype information from
    \code{reshape.gt.data}.}
  \item{x}{a data frame of genotypes from \code{fetch.gt.data}.}
  \item{by}{an expression used to define and label panels within the
    plot, evaluated in \code{data}.}
  \item{rescale}{logical: indicates if each data in each panel should
    be scaled to fill the frame.}
  \item{bounds}{contours at which to draw ellipsoid boundaries.}
  \item{min.points}{the minimum number of points for which to compute
    ellipsoids.}
  \item{between, scales, xlab, ylab, par.settings, \dots}{
    arguments passed to \code{\link[lattice:xyplot]{xyplot}}.}
}
\details{
  \code{gt.cluster.plot} is a wrapper around
  \code{\link[lattice:xyplot]{xyplot}} to facilitate plotting raw data
  underlying genotype calls.  A custom panel function draws minimal-area
  ellipsoid contours to capture proportions of genotype cluster density
  given by \code{bounds}.
  
  \code{xyplot.gt.data} is a method for packed genotype data.
}
\value{
  An object of class \code{"trellis"}.
}
\seealso{
  \code{\link{panel.cluster}},
  \code{\link[lattice:xyplot]{xyplot}},
  \code{\link[cluster:predict.ellipsoid]{ellipsoidPoints}},
  \code{\link{fetch.gt.data}},
  \code{\link{reshape.gt.data}},
  \code{\link{adjust.gt.calls}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2',raw.data=TRUE)
head(reshape.gt.data(gt))
d <- reshape.gt.data(gt[seq(32,232,40),], na.codes='n')
gt.cluster.plot(d)
gt.cluster.plot(d, scales=list(log=TRUE))
xyplot(gt[seq(32,232,40),])
}
\keyword{hplot}
