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
\name{ls.assay.group}
\alias{ls.assay.group}
\title{List Assay Groups}
\description{
  This returns a list of assay groups for the specified genotyping
  platform.
}
\usage{
ls.assay.group(platform.name, show.ids=FALSE)
}
\arguments{
  \item{platform.name}{a genotyping platform name.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per assay group, and 4 or 5 columns:
  \item{assay.group.id}{if \code{show.ids} is set: a unique integer key
    for this assay group.}
  \item{assay.group}{a short, unique identifier for the assay group.}
  \item{description}{a free-text description of the assay group.}
  \item{created.by}{the user name that created the assay group.}
  \item{created.dt}{the creation date of the assay group.}
}
\seealso{
  \code{\link{mk.assay.group}}.
}
\examples{
gt.demo.check()
ls.assay.group('Demo_Set_1')
ls.assay.group('Demo_Set_1', show.ids=TRUE)
}
\keyword{database}
