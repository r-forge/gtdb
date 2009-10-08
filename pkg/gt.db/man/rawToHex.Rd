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
\name{hexToRaw}
\alias{hexToRaw}
\alias{rawToHex}
\title{Convert between Raw Vectors and Hex Strings}
\description{
  Convert between raw vectors and hex strings.
}
\usage{
hexToRaw(hex)
rawToHex(raw)
}
\arguments{
  \item{hex}{a string of hex digits.}
  \item{raw}{a raw vector.}
}
\value{
  For \code{hexToRaw}, a raw vector constructed by converting each
  consecutive pair of hex digits to one byte.  For \code{rawToHex}, a
  character string formed by converting each byte of the input vector
  to its two-digit hex representation.
}
\seealso{
  \code{\link{charToRaw}}, \code{\link{rawToChar}}.
}
\examples{
hexToRaw('1a3b1a3b')
rawToHex(as.raw(seq(0,100,10)))
}
\keyword{classes}
