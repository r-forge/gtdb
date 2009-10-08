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
\name{as.mask}
\alias{as.mask}
\alias{un.mask}
\title{Convert To/From Character Masks}
\description{
  Convert between logical vectors, and strings of T/F characters.
}
\usage{
as.mask(x)
un.mask(s)
}
\arguments{
  \item{x}{a logical vector.}
  \item{s}{a string of T/F characters.}
}
\value{
  For \code{as.mask}, a string of T/F characters with one character
  per element of \code{x}.

  For \code{un.mask}, a logical vector with one element per character
  in \code{s}
}
\seealso{
  \code{\link{mask.str}}.
}
\examples{
as.mask(c(TRUE,FALSE,TRUE,FALSE,TRUE))
un.mask('TFFTTF')
}
\keyword{character}
