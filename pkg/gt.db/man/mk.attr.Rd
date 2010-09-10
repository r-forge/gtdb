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
\name{mk.attr}
\alias{mk.attr}
\alias{rm.attr}
\alias{ls.attr}
\title{Create, Remove, and List Attribute Definitions}
\description{
  \code{mk.attr} and \code{rm.attr} create or remove attribute
  definitions associated with various types of genotype data objects,
  such as samples and subjects, and \code{ls.attr} lists attribute
  definitions.
}
\usage{
mk.attr(target, parent.name, attr.name, datatype,
        levels, description, is.hidden=FALSE)
rm.attr(target, parent.name, attr.name)
ls.attr(target, parent.name, show.all=FALSE, show.ids=FALSE)
}
\arguments{
  \item{target}{the type of attribute to be created: i.e.,
    \code{'subject'} or \code{'sample'}.}
  \item{parent.name}{the name of the parent container for the new
    attribute, i.e. a project or dataset name.}
  \item{attr.name}{a short identifier for the new attribute.}
  \item{datatype}{the datatype for attribute values: \code{'number'},
    \code{'string'}, \code{'boolean'}, or \code{'factor'}.}
  \item{levels}{for factors, a string describing the factor levels,
    formed by concatenating a comma-separated list of single quoted
    strings.}
  \item{description}{a free text description of the attribute.}
  \item{is.hidden}{logical: indicates if the attribute is hidden.}
  \item{show.all}{logical: indicates if hidden attributes should be
    included in the output.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\details{
  These functions are usually not called directly; instead, there are
  helper functions for each type of attribute
  (i.e. \code{\link{mk.subject.attr}}, \code{\link{rm.subject.attr}},
  \code{\link{ls.subject.attr}} that are more convenient.

  Multiple attributes can be created or removed with single calls to
  \code{mk.attr} and \code{rm.attr}.  The \code{attr.name},
  \code{datatype}, \code{levels}, and \code{description} arguments can
  all be vectors and the usual recycling rules apply.
}
\value{
  For \code{mk.attr} and \code{rm.attr}, if successful, the number of
  rows inserted or deleted from the target attribute table.

  For \code{ls.attr}, a data frame with 7 columns:
  \item{attr.name}{the attribute name.}
  \item{datatype}{one of \code{'number'}, \code{'string'},
    \code{'boolean'}, or \code{'string'}.}
  \item{levels}{a single-quoted comma-separated list of factor levels.}
  \item{description}{a free text description.}
  \item{is.hidden}{logical: indicates if the attribute is hidden.}
  \item{created.by}{the user name that created the attribute.}
  \item{created.dt}{the creation date of the attribute.}
}
\seealso{
  \code{\link{mk.subject.attr}}, \code{\link{mk.sample.attr}}.
}
\examples{
\dontrun{
mk.attr('subject', 'Demo', 'stuff1', 'number',
        description='Something numeric')
mk.attr('subject', 'Demo', 'stuff2', 'factor', "'a','b','c'",
        description='A factor')
ls.attr('subject', 'Demo')
rm.attr('subject', 'Demo', c('stuff1','stuff2'))
}
}
\keyword{database}
