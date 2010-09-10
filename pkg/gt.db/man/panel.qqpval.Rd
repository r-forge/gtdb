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
\name{panel.qqpval}
\alias{prepanel.qqpval}
\alias{panel.qqpval}
\title{Quantile-Quantile Plots for P Values: Panel Functions}
\description{
  Panel functions used by \code{\link{qqpval}} for generating
  Quantile-Quantile plots of log transformed P values.
}
\usage{
prepanel.qqpval(x, n, groups=NULL, subscripts, ...)
panel.qqpval(x, n, max.pts, groups=NULL, ...)
}
\arguments{
  \item{x}{a vector of P values.}
  \item{n}{the total number of tests drawn from.}
  \item{max.pts}{the maximum number of discrete quantiles to plot.}
  \item{groups, subscripts}{see \code{\link[lattice:xyplot]{xyplot}}.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:prepanel.default]{prepanel.default.xyplot}} or
    \code{\link[lattice:panel.xyplot]{panel.qq}}.}
}
\seealso{
  \code{\link{qqpval}}, 
  \code{\link[lattice:prepanel.default]{prepanel.default.xyplot}},
  \code{\link[lattice:panel.xyplot]{panel.qq}}.
}
\keyword{hplot}
