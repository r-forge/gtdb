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
\name{unpack.gt.matrix}
\alias{unpack.gt.matrix}
\title{Convert Packed Genotype Strings to a Genotype Matrix}
\description{
  This takes a data frame of genotypes with one row per assay, and
  returns a new data frame with one column per assay and one row per
  sample.
}
\usage{
unpack.gt.matrix(gt.data, names=gt.data$assay.name, ...)
}
\arguments{
  \item{gt.data}{a data frame of genotypes from \code{fetch.gt.data}.}
  \item{names}{column names for the resulting data frame.}
  \item{\dots}{additional arguments passed to \code{\link{gt.split}}.}
}
\value{
  A data frame with one column per row in \code{gt.data}, and one row
  per sample.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{gt.split}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_1')
head(unpack.gt.matrix(gt[1:5,]))
head(unpack.gt.matrix(gt[1:5,], convert='char'))
}
\keyword{manip}
