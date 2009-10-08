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
\name{mask.gt.data}
\alias{mask.gt.data}
\title{Mask Sample Genotypes}
\description{
  Given a data frame of genotype information, this masks out a subset
  of samples based on a logical vector.
}
\usage{
mask.gt.data(gt.data, sample.mask, repack=FALSE)
}
\arguments{
  \item{gt.data}{a data frame of genotype information.}
  \item{sample.mask}{either a logical vector, or a string of T/F
    values, indicating which samples should be kept.}
  \item{repack}{logical: indicates whether masked-out samples should
    be entirely removed from the result, or left in place.}
}
\details{
  If \code{repack} is \code{FALSE}, then genotypes of masked-out
  samples are coded as missing.  If \code{repack} is \code{TRUE}, then
  masked-out samples are deleted from the result.
}
\value{
  An updated version of \code{gt.data} with genotypes for masked-out
  samples either marked as missing or deleted.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{summary.gt.data}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2')
pt <- fetch.pt.data('Demo_2')
head(summary.gt.data(gt))
gm <- mask.gt.data(gt, pt$gender=='M')
gf <- mask.gt.data(gt, pt$gender=='F')
head(summary.gt.data(gm))
head(summary.gt.data(gf))
}
\keyword{manip}
