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
\name{gt.demo.check}
\alias{gt.demo.check}
\title{Check for Presence of GT.DB Demo Datasets}
\description{
  This is used in many of the GT.DB examples to check for the presence
  of the demo datasets.
}
\usage{
gt.demo.check()
}
\details{
  If the demo datasets are unavailable, then the example code is
  terminated with an error message.  If invoked non-interactively (say,
  from \code{R CMD CHECK}), then a temporary SQLite database is created
  for the demo datasets.
}
\seealso{
  \code{\link{init.gt.db}}, \code{\link{demo}(setup.gt.demo)}.
}
\examples{
gt.demo.check()
}
\keyword{database}
