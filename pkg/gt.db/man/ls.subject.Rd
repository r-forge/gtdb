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
\name{ls.subject}
\alias{ls.subject}
\title{List Subjects in a Genotyping Project}
\description{
  This returns a list of subjects defined in the specified genotyping
  project.
}
\usage{
ls.subject(project.name, show.ids=FALSE)
}
\arguments{
  \item{project.name}{a unique project identifier.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per dataset, and 1 or 2 columns:
  \item{subject.id}{if \code{show.ids} is set: the unique integer key
    for this subject.}
  \item{subject.name}{a short, unique identifier for the subject.}
}
\seealso{
  \code{\link{ls.project}}, \code{\link{mk.subject}}.
}
\examples{
gt.demo.check()
head(ls.subject('Demo'))
head(ls.subject('Demo', show.ids=TRUE))
}
\keyword{database}
