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
\name{mask.str}
\alias{mask.str}
\title{Mask Character Strings}
\description{
  Mask selected positions from elements of a character vector.
}
\usage{
mask.str(str, mask, ch='x')
}
\arguments{
  \item{str}{a character vector to be masked.}
  \item{mask}{either a logical vector or character mask, with a 1:1
    correspondence to character positions in elements of \code{str}.}
  \item{ch}{a character that will replace masked positions in
    \code{str}.}
}
\value{
  A character vector with the same form as \code{str}, with positions
  specified as \code{FALSE} in \code{mask} either replaced by the
  first character in \code{ch}, or squeezed out entirely, if \code{ch}
  is an empty string.
}
\seealso{
  \code{\link{as.mask}}.
}
\examples{
mask.str('12345678', 'TFTFTFTF', '_')
mask.str('12345678', 'TFTFTFTF', '')
}
\keyword{character}
