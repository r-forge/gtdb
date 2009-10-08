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
\name{nsubstr}
\alias{nsubstr}
\title{Count Substring Instances}
\description{
  Counts exact instances of a target string in each element of a
  character vector.
}
\usage{
nsubstr(a,b)
}
\arguments{
  \item{a}{a character vector.}
  \item{b}{a character string with positive length.}
}
\value{
  An integer vector consisting of the number of instances of \code{b}
  in each element of \code{a}.
}
\examples{
nsubstr(c('aabbcc',NA,'aabbbbccbb'), 'bb')
}
\keyword{character}
