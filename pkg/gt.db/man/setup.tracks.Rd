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
\name{setup.tracks}
\alias{setup.tracks}
\title{Initialize Track Plot of a Genomic Interval}
\description{
  Prepare the current graphical device for drawing one or more tracks of
  information in genomic coordinates.
}
\usage{
setup.tracks(xlim, name='region', newpage=TRUE, ...)
}
\arguments{
  \item{xlim}{a vector of length 2 giving the range of genomic positions
    for the plot.}
  \item{name}{a name for the viewport defining the plotting area.}
  \item{newpage}{logical: indicates whether to erase the current
    graphical device.}
  \item{\dots}{additional arguments passed to
    \code{\link[grid:viewport]{viewport}}.}
}
\details{
  This function sets up the coordinate system to be used for drawing
  tracks.  It needs to be called before creating any tracks, because
  track construction may need to access information about the extent of
  the drawing area.
}
\seealso{
  \code{\link[grid:viewport]{viewport}}, \code{\link{draw.tracks}}.
}
\keyword{hplot}
