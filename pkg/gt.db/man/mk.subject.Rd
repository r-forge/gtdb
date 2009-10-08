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
\name{mk.subject}
\alias{mk.subject}
\alias{rm.subject}
\title{Create or Remove Project Subjects}
\description{
  \code{mk.subject} and \code{rm.subject} insert or remove subjects
  from a genotyping project.
}
\usage{
mk.subject(project.name, data)
rm.subject(project.name, subject.name)
}
\arguments{
  \item{project.name}{a short unique identifier for the project.}
  \item{data}{a data frame with a \code{subject.name} column.}
  \item{subject.name}{a vector of subject identifiers.}
}
\value{
If successful, the number of rows inserted or deleted from the subject
table.
}
\seealso{
  \code{\link{ls.subject}}, \code{\link{mk.sample}}.
}
\examples{\dontrun{
mk.subject('Demo_1', data.frame(subject.name='NA_12345'))
rm.subject('Demo_1', 'NA_12345')
}}
\keyword{database}
