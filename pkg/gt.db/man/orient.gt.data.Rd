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
\name{orient.gt.data}
\alias{orient.gt.data}
\title{Flip Assay Strands and/or Swap Alleles}
\description{
  Given a data frame of genotype information, this updates genomic
  orientations and/or swaps allele coding for individual assays.
}
\usage{
orient.gt.data(gt.data, flip=FALSE, swap=FALSE)
}
\arguments{
  \item{gt.data}{a data frame of genotype information.}
  \item{flip}{logical: indicates if strands should be flipped.}
  \item{swap}{logical: indicates if alleles should be swapped.}
}
\details{
  Both \code{flip} and \code{swap} are evaluated for each assay in
  \code{gt.data}, so these can either be single logical values or
  vectors with one value per row in \code{gt.data}.
}
\value{
  An updated version of \code{gt.data} with the specified adjustments
  applied.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{match.gt.data}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2')

# orient to dbSNP strand
gt <- orient.gt.data(gt, flip=(gt$dbsnp.orient == '-'))
with(gt, table(strand,dbsnp.orient))

# orient to '+' strand
gt <- orient.gt.data(gt, flip=(gt$strand=='-'))
with(gt, table(strand,dbsnp.orient))

# swap to major/minor allele order
swap <- (summary.gt.data(gt)$freq.a < 0.5)
table(swap)
gt <- orient.gt.data(gt, swap=swap)
table(summary.gt.data(gt)$freq.a < 0.5)

# sort alleles
head(gt$alleles)
swap <- (sapply(strsplit(gt$alleles,'/'), order)[1,] == 2)
gt <- orient.gt.data(gt, swap=swap)
head(gt$alleles)
}
\keyword{manip}
