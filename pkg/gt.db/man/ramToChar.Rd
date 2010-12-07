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
\name{charToRam}
\alias{charToRam}
\alias{ramToChar}
\title{Convert between Raw Matrices and Character Vectors}
\description{
  Convert between raw matrices and character vectors.
}
\usage{
charToRam(str)
ramToChar(raw)
}
\arguments{
  \item{str}{a character vector.}
  \item{raw}{a raw matrix.}
}
\value{
  For \code{charToRam}, a raw matrix with one column per string in
  \code{char}.
  
  For \code{ramToChar}, a character vector with one element per column
  in \code{raw}.
}
\seealso{
  \code{\link{hexToRam}}, \code{\link{ramToHex}},
  \code{\link{charToRaw}}, \code{\link{rawToChar}}.
}
\examples{
charToRam(c('1234','5678','abcd'))
m <- matrix(charToRaw('abcd'),4,3)
print(m)
ramToChar(m)
}
\keyword{classes}
