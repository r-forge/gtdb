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
\name{qqpval}
\alias{qqpval}
\title{Quantile-Quantile Plots for P Values}
\description{
  Quantile-Quantile plots of log transformed P values versus their
  expected uniform distribution, with support for very large numbers
  of tests.
}
\usage{
qqpval(x, ..., n=NA, max.pts=1000, prepanel=prepanel.qqpval,
       panel=panel.qqpval, xlab='theoretical quantiles')
}
\arguments{
  \item{x}{either a vector of P values, or a formula to be evaluated
    in the context of a \code{data} argument.}
  \item{\dots}{additional arguments passed to
    \code{\link[lattice:qqmath]{qqmath}}.}
  \item{n}{the total number of tests drawn from.}
  \item{max.pts}{the maximum number of discrete quantiles to plot.}
  \item{prepanel, panel}{as in \code{\link[lattice:qqmath]{qqmath}}.}
  \item{xlab}{X axis label.}
}
\details{
  Conventional Q-Q plots for very large numbers of tests are slow and
  most of the \dQuote{ink} is spent plotting uninteresting large P
  values.  This function avoids this problem by limiting the number of
  quantiles plotted.  If the number of supplied P values is larger than
  \code{max.pts}, then a threshold is chosen to switch between plotting
  each observed value, and plotting evenly spaced quantiles.  The
  threshold is optimized to minimize the quantile spacing in the
  equal-spaced region.

  A valid (truncated) Q-Q plot can also be generated from the smallest P
  values selected from a larger set of tests, by specifying \code{n},
  the total number of tests the supplied values were selected from.

  A reference line is drawn using settings specified by
  \code{trellis.par.get('reference.line')}.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{qqthin}}, \code{\link[lattice:qqmath]{qqmath}}.
}
\examples{
pval <- runif(250000)
qqpval(pval)
qqpval(pval, max.pts=100)
qqpval(pval[pval<0.001], n=250000, max.pts=100)

p <- runif(3000)
s <- replicate(400, sort(runif(3000)))
q <- apply(s, 1, quantile, c(0.025,0.975))
rc <- '#b0b0b0'
qqpval(~q[1,]+q[2,]+p, type=c('l','p'), col=c(rc,rc,2),
       lty=c(2,2,0), pch=c(32,32,1), max.pts=100,
       par.settings=list(reference.line=list(col=rc)))
}
\keyword{hplot}
