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
\name{panel.cluster}
\alias{panel.cluster}
\title{Panel Function for Drawing Elliptical Cluster Boundaries}
\description{
  Panel function for plotting two-dimensional clusters of points with
  elliptical boundary regions.
}
\usage{
panel.cluster(x, y, group.number=1, bounds=c(), min.points=4, ...)
}
\arguments{
  \item{x, y}{as in \code{\link[lattice:panel.xyplot]{panel.xyplot}}.}
  \item{group.number}{as in
    \code{\link[lattice:panel.superpose]{panel.superpose}}.}
  \item{bounds}{contours at which to draw ellipsoid boundaries.}
  \item{min.points}{the minimum number of points for which to compute
    ellipsoids.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:panel.xyplot]{panel.xyplot}}.}
}
\details{
  Points are first plotted as in
  \code{\link[lattice:panel.xyplot]{panel.xyplot}}.  Then
  elliptical boundaries are drawn using settings from
  \code{superpose.line}, using \code{trellis.par.get}.
}
\seealso{
  \code{\link{gt.cluster.plot}},
  \code{\link[lattice:panel.xyplot]{panel.xyplot}},
  \code{\link[lattice:panel.superpose]{panel.superpose}},
  \code{\link[cluster:predict.ellipsoid]{ellipsoidPoints}}.
}
\examples{
g <- rep(0:1,100)
x <- rnorm(200)
y <- g*6 + 2*g*x + rnorm(100)
xyplot(y~x, groups=g, panel=panel.superpose,
       panel.groups=panel.cluster, bounds=c(0.5,0.8,0.95))
}
\keyword{hplot}
