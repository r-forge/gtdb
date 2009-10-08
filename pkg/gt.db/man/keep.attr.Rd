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
\name{keep.attr}
\alias{keep.attr}
\alias{kept.attr}
\alias{[.keep.attr}
\alias{[<-.keep.attr}
\title{Keep User Attributes}
\description{
  This function creates objects with the property that they more
  systematically preserve user attributes when they are indexed.
}
\usage{
keep.attr(.Data, ..., .Attr=NULL)
kept.attr(x)
\method{[}{keep.attr}(x, \dots)
\method{[}{keep.attr}(x, \dots) <- value
}
\arguments{
  \item{.Data}{an R object for which attributes are to be kept.}
  \item{...}{attributes in key=value form to be attached to the object.}
  \item{.Attr}{additional attributes to attach, collected into a list.}
  \item{x}{an object of class \code{'keep.attr'}.}
  \item{value}{a suitable replacement value.}
}
\value{
  For \code{keep.attr}, the original object with additional class
  \code{'keep.attr'}, with additional attributes attached.  Indexing
  on this object will propagate all user specified attributes, so long
  as the indexing operation returns an object of the same class.
  Thus, for a data frame, attributes will be preserved for operations
  that return a new data frame.

  For \code{kept.attr}, a list of all the user defined attributes
  associated with the object.
}
\seealso{\code{\link{attributes}}.}
\examples{
d <- data.frame(a=c(1,2,3), b=c(4,5,6), c=c(7,8,9))
d <- structure(d, xyz='something')
e <- keep.attr(d, abc='more')
d[3] <- 5
e[3] <- 5
str(d)
str(d[1])
str(d[,-1])
str(d[,1])
str(e)
str(e[1])
str(e[,-1])
str(e[,1])
}
\keyword{classes}
