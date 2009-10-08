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
\name{mk.project}
\alias{mk.project}
\alias{rm.project}
\title{Create or Remove a Genotyping Project}
\description{
  These functions insert or remove entries from the project table.
}
\usage{
mk.project(project.name, description, is.hidden=FALSE)
rm.project(project.name)
}
\arguments{
  \item{project.name}{a short unique identifier for the project.}
  \item{description}{a free-text description of the project.}
  \item{is.hidden}{logical: indicates if the project is hidden.}
}
\value{
If successful, the number of rows inserted or removed from the project
table (i.e., 1). 
}
\seealso{
  \code{\link{ls.project}}, \code{\link{set.hidden}}.
}
\examples{
\dontrun{
mk.project('Demo_9', 'Demo Genome-Wide Association Study')
rm.project('Demo_9')
}
}
\keyword{database}
