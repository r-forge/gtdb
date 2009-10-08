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
\name{ch.table}
\alias{ch.table}
\title{Character Based Contingency Table}
\description{
  Build contingency tables of character counts.
}
\usage{
ch.table(s1, s2, chars)
}
\arguments{
  \item{s1}{a vector of character strings.}
  \item{s2}{an optional vector of character strings, with lengths
    matching strings in \code{s1}.}
  \item{chars}{a vector of individual characters to tabulate.}
}
\value{
  If \code{s2} is missing, a matrix of character counts with one row
  per element of \code{s1} and one column per element of \code{chars}.

  If \code{s2} is present, a three dimensional array where the first
  dimension indicates a position within \code{s1} and \code{s2} (with
  recycling if necessary, and the second and third dimensions form
  contingency tables of corresponding characters in these strings.
}
\examples{
s1 <- c('AACAGCTACAGT','TTGTCGATGTCA')
s2 <- 'AACAGCTACAGT'
ch.table(s1, chars=c('A','C','G','T'))
x <- ch.table(s1, s2, chars=c('A','C','G','T'))
x[1,,]
x[2,,]
}
\keyword{character}
