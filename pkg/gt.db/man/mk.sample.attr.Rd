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
\name{mk.sample.attr}
\alias{mk.sample.attr}
\alias{rm.sample.attr}
\alias{ls.sample.attr}
\title{Create, Remove, or List Sample Attribute Definitions}
\description{
  \code{mk.sample.attr} and \code{rm.sample.attr} create or remove
  sample attributes associated with a particular genotype dataset, and
  \code{ls.sample.attr} lists sample attributes for this dataset.
}
\usage{
mk.sample.attr(dataset.name, data, description=names(data))
rm.sample.attr(dataset.name, attr.name)
ls.sample.attr(dataset.name)
}
\arguments{
  \item{dataset.name}{the unique identifier for a dataset.}
  \item{data}{a data frame with one or more columns to be used to
    model the new attributes.  See details.}
  \item{description}{a vector of free-text descriptions with one
    element per column in \code{data}.}
  \item{attr.name}{a vector of attribute names to be removed.}
}
\details{
  These functions are wrappers around \code{\link{mk.attr}},
  \code{\link{rm.attr}}, and \code{\link{ls.attr}}.  A call to
  \code{mk.sample.attr} results in creation of a new attribute for
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
sd <- fetch.sample.data('Demo_1')
sd$stuff.1 <- rep(1:3, length.out=nrow(sd))
sd$stuff.2 <- factor(sd$stuff.1, levels=1:3, labels=c('a','b','c'))
sd$stuff.3 <- as.character(sd$stuff.2)
str(sd)
mk.sample.attr('Demo_1', sd[c('stuff.1','stuff.2','stuff.3')])
ls.sample.attr('Demo_1')
rm.sample.attr('Demo_1', c('stuff.1','stuff.2','stuff.3'))
}
}
\keyword{database}
