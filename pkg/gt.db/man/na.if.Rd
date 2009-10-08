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
\name{na.if}
\alias{na.if}
\title{Conditional Conversion to Missing Values}
\description{
  \code{na.if} returns a value with the same shape as \code{v1},
  with values equal to \code{v2} set to \code{NA}.
}
\usage{
na.if(v1, v2)
}
\arguments{
  \item{v1}{an object to be tested.}
  \item{v2}{a value to be compared to elements of \code{v1}.}
}
\details{
  The result is equivalent to \code{ifelse(v1==v2, NA, v1)}.
}
\value{
  An object of the same shape as \code{v1}, where values equal to
  \code{v2} are set to \code{NA}.
}
\references{
  Inspired by Oracle's \code{nullif()} function.
}
\seealso{
  \code{\link{if.na}}, \code{\link{ifelse}}.
}
\examples{
x <- c(1:3,1:3)
na.if(x, 3)
}
\keyword{logic}
\keyword{NA}
