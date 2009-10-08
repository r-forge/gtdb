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
\name{ls.mapping}
\alias{ls.mapping}
\title{List Assay Mapping Sets}
\description{
  This returns a list of mapping sets for the specified genotyping
  platform.
}
\usage{
ls.mapping(platform.name, show.all=FALSE, show.ids=FALSE)
}
\arguments{
  \item{platform.name}{a genotyping platform name.}
  \item{show.all}{logical: indicates if hidden mapping sets should be
    included in the output.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per mapping, and 4 or 5 columns:
  \item{mapping.id}{if \code{show.ids} is set: a unique integer key
    for this mapping.}
  \item{mapping}{a short, unique identifier for the mapping.}
  \item{description}{a free-text description of the mapping.}
  \item{created.by}{the user name that created the mapping.}
  \item{created.dt}{the creation date of the mapping.}
}
\seealso{
  \code{\link{mk.mapping}}.
}
\examples{
gt.demo.check()
ls.mapping('Demo_Set_1')
ls.mapping('Demo_Set_1', show.ids=TRUE)
}
\keyword{database}
