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
\name{gplot.prcomp}
\alias{gplot.prcomp}
\title{Genome-Wide Level Plots of SNP Loadings from PCA}
\description{
  Generates level plots of genome-wide SNP loadings from a principal
  components analysis.
}
\usage{
gplot.prcomp(x, col=1:5, aggr.fn=function(x) sum(x^2),
             rescale=TRUE, xlab=NULL, ylab=NULL, ...)
}
\arguments{
  \item{x}{a data structure returned by \code{\link{snp.loadings}}.}
  \item{col}{a vector of component numbers to plot.}
  \item{aggr.fn}{an aggregating function to apply to data within each bin.}
  \item{rescale}{logical: indicates if results should be scaled based
    on the genome-wide average bin value.}
  \item{xlab, ylab}{axis labels passed to \code{\link{gplot}}.}
  \item{\dots}{additional arguments passed to \code{\link{gplot}}.}
}
\details{
  The output of \code{gplot.prcomp} is a vertically stacked series of
  genome-wide plots of SNP loadings for a specified set of principal
  components.  The default aggregating function results in a plot of
  variance binned by genomic interval.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{prcomp.gt.data}}, \code{\link{prcomp.gt.dataset}},
  \code{\link{snp.loadings}}, \code{\link{gplot}}.
}
\examples{
gt.demo.check()
p <- prcomp(gt.dataset('Demo_1', gt.filter=(ploidy=='A')))
s <- snp.loadings(p, gt.dataset('Demo_1', gt.filter=TRUE))
par <- list(fontsize=list(text=8),
            background=list(alpha=0,col=gray(0.9)))
gplot.prcomp(s, col=1:3, par.settings=par)
}
\keyword{hplot}
