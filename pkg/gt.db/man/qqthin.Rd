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
\name{qqthin}
\alias{qqthin}
\title{Sparse Normal Quantile-Quantile Plots}
\description{
  Quantile-Quantile plots of a sample versus the normal distribution,
  with support for very large numbers of observations.
}
\usage{
qqthin(x, ..., max.pts=1000, panel=panel.qqthin)
}
\arguments{
  \item{x}{either a numeric vector, or a formula to be evaluated
    in the context of a \code{data} argument.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:qqmath]{qqmath}}.}
  \item{max.pts}{the maximum number of discrete quantiles to plot.}
  \item{panel}{as in \code{\link[lattice:qqmath]{qqmath}}.}
}
\details{
  Conventional Q-Q plots for very large numbers of observations are slow
  and most of the \dQuote{ink} is spent plotting uninteresting values.
  This function avoids this problem by limiting the number of distinct
  quantiles plotted.  If the number of values is larger than
  \code{max.pts}, then a threshold is chosen to switch between plotting
  each extreme value, and plotting evenly spaced quantiles.  The
  threshold is optimized to minimize the quantile spacing in the
  equal-spaced region.

  A reference line is drawn using settings specified by
  \code{trellis.par.get('reference.line')}.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{qqpval}}, \code{\link[lattice:qqmath]{qqmath}}.
}
\examples{
v <- rnorm(100000)
qqthin(v)
qqthin(v, max.pts=100)

v <- rnorm(3000)
s <- replicate(400, sort(rnorm(3000)))
q <- apply(s, 1, quantile, c(0.025,0.975))
rc <- '#b0b0b0'
qqthin(~q[1,]+q[2,]+v, type=c('l','p'), col=c(rc,rc,2),
       lty=c(2,2,0), pch=c(32,32,1), max.pts=100,
       par.settings=list(reference.line=list(col=rc)))
}
\keyword{hplot}
