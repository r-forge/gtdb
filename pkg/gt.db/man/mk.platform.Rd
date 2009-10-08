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
\name{mk.platform}
\alias{mk.platform}
\alias{rm.platform}
\title{Create or Remove a Genotyping Platform}
\description{
  These functions insert or remove entries from the genotyping
  platform table. 
}
\usage{
mk.platform(platform.name, description)
rm.platform(platform.name)
}
\arguments{
  \item{platform.name}{a short unique identifier for the platform.}
  \item{description}{a free-text description of the platform.}
}
\value{
If successful, the number of rows inserted in or removed from the
platform table (i.e., 1). 
}
\seealso{
  \code{\link{ls.platform}}, \code{\link{mk.assay.group}}.
}
\examples{\dontrun{
mk.platform('GT_800K', 'My 800K Genotyping Platform')
rm.platform('GT_800K')
}}
\keyword{database}
