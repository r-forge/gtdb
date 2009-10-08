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
\name{nw.orient.assay}
\alias{nw.orient.assay}
\title{Orient Assay Sequences}
\description{
  Determines the relative orientation of pairs of genotyping assays
  by sequence alignment.
}
\usage{
nw.orient.assay(gt1, gt2, delta=1)
}
\arguments{
  \item{gt1, gt2}{assay sequence(s) with a variant position denoted by
    an underscore (\code{'_'}).}
  \item{delta}{the minimum difference in scores to establish orientation.}
}
\details{
  Corresponding elements of \code{gt1} and \code{gt2} are split into
  left and right flanks around the variant position, and then aligned
  in forward and reverse-complement orientations.  The best-scoring
  orientation is reported if its score exceeds the score of the
  opposing orientation by at least \code{delta}.
}
\value{
  A vector of \code{'+'}, \code{'-'}, or \code{NA} values.
}
\seealso{
  \code{\link{nw.align}}.
}
\examples{
a1 <- 'ACAAC_ATGCT'
a2 <-  'GCAT_GTTG'
nw.orient.assay(a1,a2)
}
\keyword{manip}
