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
\name{if.na}
\alias{if.na}
\title{Conditional Element Selection for Missing Values}
\description{
  \code{if.na} returns a value with the same shape as \code{val},
  populated from either \code{yes} or \code{no} depending on whether
  the element of \code{is.na(val)} is \code{TRUE} or \code{FALSE}.
}
\usage{
if.na(val, yes, no=val)
}
\arguments{
  \item{val}{an object to be tested for \code{NA} values.}
  \item{yes}{return values for \code{NA} elements of \code{val}.}
  \item{no}{return values for non-\code{NA} elements of \code{val}.}
}
\details{
  The result is equivalent to \code{ifelse(is.na(val), yes, no)}.
}
\value{
  A vector of the same length as \code{val}, where missing values in
  are taken from \code{yes}, and non-missing values are taken from
  \code{no} (or \code{val} itself, if \code{no} is missing). 
}
\references{
  Inspired by Oracle's \code{nvl()} and \code{nvl2()} functions.
}
\seealso{
  \code{\link{is.na}}, \code{\link{ifelse}}.
}
\examples{
x <- c(1:3,NA,NA)
if.na(x, 4)
if.na(x, 0, x+1)
}
\keyword{logic}
\keyword{NA}
