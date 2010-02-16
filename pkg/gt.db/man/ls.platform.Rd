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
\name{ls.platform}
\alias{ls.platform}
\title{List Genotyping Platforms}
\description{
  This returns a list of genotyping platforms defined in the current
  database. 
}
\usage{
ls.platform(platform.name='\%', show.ids=FALSE)
}
\arguments{
  \item{platform.name}{an SQL \code{LIKE} expression for matching
    platform names.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per platform, and 5 or 6 columns:
  \item{platform.id}{if \code{show.ids} is set: a unique integer key
    for this platform.}
  \item{platform.name}{a short, unique identifier for the platform.}
  \item{description}{a free-text description of the platform.}
  \item{datasets}{a count of the number of datasets using this platform.}
  \item{created.by}{the user name that created the platform.}
  \item{created.dt}{the creation date of the platform.}
}
\seealso{
  \code{\link{mk.platform}}.
}
\examples{
gt.demo.check()
ls.platform()
ls.platform(show.ids=TRUE)
}
\keyword{database}
