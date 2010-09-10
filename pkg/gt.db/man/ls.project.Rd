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
\name{ls.project}
\alias{ls.project}
\title{List Genotyping Projects}
\description{
  This returns a list of genotyping projects defined in the current
  database. 
}
\usage{
ls.project(project.name='\%', show.all=FALSE, show.ids=FALSE)
}
\arguments{
  \item{project.name}{an SQL \code{LIKE} expression for matching
    project names.}
  \item{show.all}{logical: indicates if hidden projects should be
    included in the output.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per project, and 6 or 7 columns:
  \item{project.id}{if \code{show.ids} is set: a unique integer key
    for this project.}
  \item{project.name}{a short, unique identifier for the project.}
  \item{description}{a free-text description of the project.}
  \item{datasets}{a count of the number of datasets under this project.}
  \item{is.hidden}{logical: indicates if the project is hidden.}
  \item{created.by}{the user name that created the project.}
  \item{created.dt}{the creation date of the project.}
}
\seealso{
  \code{\link{mk.project}}.
}
\examples{
gt.demo.check()
ls.project()
ls.project(show.ids=TRUE)
}
\keyword{database}
