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
\name{panel.qqthin}
\alias{panel.qqthin}
\title{Sparse Normal Quantile-Quantile Plots: Panel Function}
\description{
  Panel function used by \code{\link{qqthin}} for generating
  quantile-quantile plots of very large numbers of observations.
}
\usage{
panel.qqthin(x, max.pts, groups=NULL, ...)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{max.pts}{the maximum number of discrete quantiles to plot.}
  \item{groups}{see \code{\link[lattice:xyplot]{xyplot}}.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:panel.xyplot]{panel.qqmath}}.}
}
\seealso{
  \code{\link{qqthin}}, 
  \code{\link[lattice:panel.xyplot]{panel.qqmath}}.
}
\keyword{hplot}
