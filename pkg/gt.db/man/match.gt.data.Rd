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
\name{match.gt.data}
\alias{match.gt.data}
\title{Identify Equivalent Genotyping Assays}
\description{
  This identifies assays from two genotype datasets that correspond to
  the same polymorphisms, based either on position or dbSNP rsIDs.
}
\usage{
match.gt.data(gt.data.1, gt.data.2, by=c('position','dbsnp.rsid'),
              all=FALSE)
}
\arguments{
  \item{gt.data.1, gt.data.2}{data frames with assay information to be
    compared.}
  \item{by}{the method to use to identify equivalent assays.}
  \item{all}{logical: indicates whether to report assays that could not
    be oriented and/or had inconsistent alleles.}
}
\details{
  This first identifies assays with corresponding positions (or dbSNP
  rsIDs) across the two datasets.  For each matching pair, the
  reported genomic orientation and alleles are checked for
  compatibility.  Putative matches are rejected if the alleles are
  inconsistent.
}
\value{
  A data frame with four columns.  The first two columns contain
  corresponding assay names from the two datasets.  The remaining two
  columns are:
  \item{is.flipped}{logical: indicates that the assays are in
    opposing genomic orientations.}
  \item{is.swapped}{logical: indicates that the A/B alleles are
    interchanged between assays.}
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{orient.gt.data}}.
}
\examples{
gt.demo.check()
g1 <- fetch.gt.data('Demo_2')
g2 <- orient.gt.data(g1, flip=(g1$dbsnp.orient == '-'))
head(match.gt.data(g1,g2))
}
\keyword{manip}
