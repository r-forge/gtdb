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
\name{set.hidden}
\alias{set.hidden}
\title{Update Hidden Status}
\description{
  This updates the \dQuote{hidden} status of a database object
  (project, dataset, platform mapping, sample attribute, etc).
}
\usage{
set.hidden(table, name, is.hidden=TRUE, ...)
}
\arguments{
  \item{table}{the type of object to be hidden (or revealed).}
  \item{name}{the name of the object to be hidden (or revealed).}
  \item{is.hidden}{logical: the new status for the object.}
  \item{\dots}{additional qualifiers to identify the target object.}
}
\value{
If successful, the number of rows updated in the corresponding table
(i.e., 1).
} 
\seealso{
  \code{\link{mk.project}}, \code{\link{mk.dataset}}.
}
\examples{
\dontrun{
set.hidden('dataset','Demo_1',is.hidden=TRUE)
ls.dataset('Demo')
set.hidden('dataset','Demo_1',is.hidden=FALSE)
ls.dataset('Demo')
set.hidden('subject_attr','plate',is.hidden=TRUE,project.name='Demo')
}
}
\keyword{database}
