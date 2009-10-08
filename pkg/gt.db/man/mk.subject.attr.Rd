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
\name{mk.subject.attr}
\alias{mk.subject.attr}
\alias{rm.subject.attr}
\alias{ls.subject.attr}
\title{Create, Remove, or List Subject Attribute Definitions}
\description{
  \code{mk.subject.attr} and \code{rm.subject.attr} create or remove
  subject attributes associated with a particular genotyping project,
  and \code{ls.subject.attr} lists subject attributes for this project.
}
\usage{
mk.subject.attr(project.name, data, description=names(data))
rm.subject.attr(project.name, attr.name)
ls.subject.attr(project.name)
}
\arguments{
  \item{project.name}{the unique identifier for a project.}
  \item{data}{a data frame with one or more columns to be used to
    model the new attributes.  See details.}
  \item{description}{a vector of free-text descriptions with one
    element per column in \code{data}.}
  \item{attr.name}{a vector of attribute names to be removed.}
}
\details{
  These functions are wrappers around \code{\link{mk.attr}} and
  \code{\link{rm.attr}}, and \code{\link{ls.attr}}.  A call to
  \code{mk.subject.attr} results in creation of a new attribute for
  each column in the \code{data} argument, based on the data type,
  factor levels, etc of that column.
}
\value{
  For \code{mk.sample.attr} and \code{rm.sample.attr}, the number of
  attribute definitions inserted or removed from the sample attribute
  table. 

  For \code{ls.sample.attr}, a data frame with 7 columns describing
  subject attributes defined for the specified project, as defined in
  \code{\link{ls.attr}}.
}
\seealso{
  \code{\link{mk.attr}}, \code{\link{rm.attr}}, \code{\link{ls.attr}}.
}
\examples{
\dontrun{
sd <- fetch.subject.data('Demo')
sd$stuff.1 <- rep(1:3, length.out=nrow(sd))
sd$stuff.2 <- factor(sd$stuff.1, levels=1:3, labels=c('a','b','c'))
sd$stuff.3 <- as.character(sd$stuff.2)
str(sd)
mk.subject.attr('Demo', sd[c('stuff.1','stuff.2','stuff.3')])
ls.subject.attr('Demo')
rm.subject.attr('Demo', c('stuff.1','stuff.2','stuff.3'))
}
}
\keyword{database}
