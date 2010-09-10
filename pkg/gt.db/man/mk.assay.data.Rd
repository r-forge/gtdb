%
% Copyright (C) 2009, Perlegen Sciences, Inc.
% Copyright (C) 2010, 23andMe, Inc.
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
\name{mk.assay.data}
\alias{mk.assay.data}
\alias{rm.assay.data}
\title{Create Genotyping Assay Data}
\description{
  \code{mk.assay.data} defines genotype data associated with already
  defined genotype assays.
}
\usage{
mk.assay.data(dataset.name, data, progress=FALSE)
}
\arguments{
  \item{dataset.name}{a short unique identifier for the genotype dataset.}
  \item{data}{a data frame with one row per assay.  See details.}
  \item{progress}{logical: indicates whether to report progress while
    the data is loaded.}
}
\details{
  A data frame of genotype data can provide up to five columns of
  information: 
  \describe{
    \item{assay.id}{the unique integer ID for this assay.}
    \item{flags}{an integer value composed of single-bit flags.}
    \item{genotype}{a string of a/h/b/n/x genotypes.}
    \item{qscore}{a hex string of packed quality scores.}
    \item{raw.data}{a hex string containing packed raw data (such as
      signal intensities or read counts).}
  }
  An \code{assay.name} column can be supplied in place of
  \code{assay.id}.

  Three raw data layouts are currently supported.  The 'signal' layout
  consists of a 16-bit unsigned little endian value for each allele.
  The 'seqread' layout consists of 8-bit unsigned counts for forward and
  reverse orientations for each allele, representing read counts from a
  sequencing experiment.  The 'chpdata' layout includes a 16-bit 8.8
  fixed point log allele ratio, a 16-bit 8.8 log signal strength, and 
  a 'forced' genotype call.
}
\value{
  If successful, the number of rows inserted into the assay data
  table.
}
\seealso{
  \code{\link{mk.assay.position}},
  \code{\link{mk.assay}},
  \code{\link{fetch.gt.data}},
  \code{\link{pack.raw.data}},
  \code{\link{unpack.raw.data}}.
}
\examples{\dontrun{
data(demo_01)
head(assay.pos.01)
mk.mapping('Demo_Set_1', 'Demo_Set_1_b36',
           'Platform 1 Mapping to NCBI Build 36')
mk.assay.position('Demo_Set_1', data=assay.pos.01)
}}
\keyword{database}
