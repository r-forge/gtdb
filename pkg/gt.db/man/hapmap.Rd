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
\name{hapmap}
\docType{data}
\alias{hapmap.subjects}
\title{Subject Data from the International HapMap Project}
\description{
  \code{hapmap.subjects} describes sample plate and panel membership,
  and parent/child relationships, for the 1301 individuals genotyped in
  the Phase II and Phase III HapMap Projects.
}
\usage{data(hapmap)}
\source{The International HapMap Project}
\references{
  The International HapMap Project. \url{http://www.hapmap.org}.
}
\seealso{
  \code{\link{demo}(setup.hapmap)}.
}
\examples{\dontrun{
# create HapMap project in GT.DB
demo('setup.hapmap')
}}
\keyword{datasets}
