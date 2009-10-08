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
\name{mk.assay.position}
\alias{mk.assay.position}
\alias{rm.assay.position}
\title{Create Genotyping Assay Positions}
\description{
  \code{mk.assay.position} defines map positions associated with
  genotyping assays.
}
\usage{
mk.assay.position(platform.name, mapping.name, data)
}
\arguments{
  \item{platform.name}{a short unique identifier for the platform.}
  \item{mapping.name}{the platform mapping to be populated.}
  \item{data}{a data frame with one row per assay.  See details.}
}
\details{
  If \code{mapping.name} is missing, it will default to the current
  (visible) mapping for the specified platform, if that is unique.
  A data frame of assay positions can provide up to seven columns of
  information: 
  \describe{
    \item{assay.name}{the assay identifier.}
    \item{scaffold}{a string identifying the sequence to which the
      assay is mapped.}
    \item{position}{a one-based position within the specified scaffold.}
    \item{strand}{either \code{"+"}, \code{"-"}, or \code{NA}.}
    \item{ploidy}{describes the expected allele count: it can take
      values \code{"A"} (autosomal), \code{"M"} (mitochondrial, or
      haploid), \code{"X"} (diploid in females, haploid in males), or
      \code{"Y"} (haploid in males, absent in females).  Pseudoautosomal
      positions are coded as \code{"A"}.}
    \item{dbsnp.rsid}{the dbSNP refSNP cluster ID for this assay.}
    \item{dbsnp.orient}{the orientation of the assay compared to the
      dbSNP cluster: either \code{"+"}, \code{"-"}, or \code{NA}.}
  }
}
\value{
  If successful, the number of rows inserted into the assay map
  position table.
}
\seealso{
  \code{\link{ls.assay.position}},
  \code{\link{mk.mapping}},
  \code{\link{mk.assay}}.
}
\examples{\dontrun{
data(demo_01)
head(assay.pos.01)
mk.mapping('Demo_Set_1', 'Demo_Set_1_b36',
           'Platform 1 Mapping to NCBI Build 36')
mk.assay.position('Demo_Set_1', data=assay.pos.01)
}}
\keyword{database}
