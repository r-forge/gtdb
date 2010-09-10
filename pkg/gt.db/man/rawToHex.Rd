%
% Copyright (C) 2009, Perlegen Sciences, Inc.
% Copyright (C) 2010, 23andMe, Inc.
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
\alias{hexToRam}
\alias{ramToHex}
\title{Convert between Raw Vectors and Hex Strings}
\description{
  Convert between raw vectors and hex strings.
}
\usage{
hexToRaw(hex)
rawToHex(raw)
hexToRam(hex)
ramToHex(raw)
}
\arguments{
  \item{hex}{a character vector composed of strings of hex digit pairs,
    all the same length.}
  \item{raw}{a raw vector or matrix.}
}
\details{
  \code{hexToRaw} and \code{RawToHex} convert between a single character
  string and a raw vector.  \code{hexToRam} and \code{ramToHex} convert
  between a character vector and a raw matrix.  Columns of the raw
  matrix become individual character strings, and vice versa.
}
\value{
  For \code{hexToRaw}, a raw vector with one byte per pair of hex digits 
  in the input string.  For \code{hexToRam}, a raw matrix with one
  column per string in \code{hex}, constructed by converting each
  consecutive pair of hex digits to one byte.
  
  For \code{rawToHex}, a character string formed by converting each byte
  of the input vector to its two-digit hex representation.  For
  \code{ramToHex}, a vector formed by converting each column of an input
  matrix to its hex representation.
}
\seealso{
  \code{\link{charToRaw}}, \code{\link{rawToChar}},
  \code{\link{charToRam}}, \code{\link{ramToChar}}.
}
\examples{
hexToRaw('1a3b1a3b')
rawToHex(as.raw(seq(0,100,10)))
hexToRam(c('1234','5678','9abc'))
}
\keyword{classes}
