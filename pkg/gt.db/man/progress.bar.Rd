%
% Copyright (C) 2009, Perlegen Sciences, Inc.
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
\name{progress.bar}
\alias{progress.bar}
\title{Console Text-Based Progress Bar}
\description{
  Displays a text-based progress bar on the R console, to indicate
  how much of a long-running calculation has been completed.  The
  elapsed time and estimated time remaining are also reported.
}
\usage{
progress.bar(done, total, width=getOption('width'))
}
\arguments{
  \item{done}{the amount of work already done.}
  \item{total}{the total amount of work to do.}
  \item{width}{the console width, in characters.}
}
\details{
  After writing the progress bar, the cursor is repositioned at the
  start of the line so that successive calls will update the bar in
  place.  The console buffer is flushed after each call.  The first
  call should have \code{done=0} to ensure proper initialization of
  internal data structures.
}
\value{
  Invisible \code{NULL}.
}
\examples{
fn <- function(x) { Sys.sleep(1); progress.bar(x, 20) }
x <- sapply(0:20, fn)
}
\keyword{utilities}
