%
% Copyright (C) 2010, 23andMe, Inc.
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
\name{draw.tracks}
\alias{draw.tracks}
\title{Draw Tracks in Genomic Coordinates}
\description{
  Draw one or more tracks of genomically structured information in the
  current graphical context.
}
\usage{
draw.tracks(tracks, scale=c('Mb','Kb','bp'), xlab)
}
\arguments{
  \item{tracks}{a list of track objects to draw in the current
    viewport, arranged from bottom to top.}
  \item{scale}{a scaling to use for annotating the X axis of the track
    plot.}
  \item{xlab}{a label for the X axis.}
}
\details{
  This function draws a series of one or more tracks of data organized
  in genomic coordinates within the current viewport, as initialized by
  \code{\link{setup.tracks}}.  Each track is drawn in its own viewport;
  if the elements of \code{tracks} are named, then these names are used
  to name the corresponding viewports.
}
\seealso{
  \code{\link{setup.tracks}}, \code{\link{manhattan.track}},
  \code{\link{gene.track}}, \code{\link{ld.track}}.
}
\examples{
gt.demo.check()
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
s  <- score.gt.data(status~genotype, pt, gt)
setup.tracks(range(gt$position), width=unit(4,'inch'),
             y=unit(1,'inch'), just='bottom')
mt <- manhattan.track(gt[rownames(s),]$position, -log10(s$pvalue))
ld <- ld.track(gt)
draw.tracks(list(ld=ld, mt=mt))
}
\keyword{hplot}
