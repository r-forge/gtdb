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
\name{revcomp}
\alias{revcomp}
\title{Reverse Complement DNA Sequences}
\description{
  Reverse complement DNA sequences, including IUPAC ambiguity codes.
}
\usage{
revcomp(x)
}
\arguments{
  \item{x}{a character vector of DNA sequences}
}
\value{
  A vector with the same shape as \code{x}, where each sequence has been
  replaced by its reverse complement, preserving case.
}
\examples{
revcomp(c('AACAGTAGA','AACCNNRacgt'))
}
\keyword{manip}
