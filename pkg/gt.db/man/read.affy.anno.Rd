%
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
\name{read.affy.anno}
\alias{read.affy.anno}
\alias{load.affy.platform}
\alias{load.affy.mapping}
\title{Import Affymetrix NetAffx Annotation Data}
\description{
  The \code{read.affy.anno} function reads the contents of an Affymetrix
  annotation CSV file into a data frame; \code{load.affy.platform} loads
  assay definitions (flanks, alleles) into the database; and
  \code{load.affy.mapping} loads assay positions and dbSNP information.
}
\usage{
read.affy.anno(file)
load.affy.platform(anno, description, progress=TRUE)
load.affy.mapping(anno, progress=TRUE)
}
\arguments{
  \item{file}{the file name of the annotation file.}
  \item{anno}{an annotation data frame from \code{read.affy.anno}.}
  \item{description}{a free-text description of the platform.}
  \item{progress}{logical: indicates whether to report progress during
    the database load.}
}
\details{
  \code{load.affy.platform} creates and populates a new platform, using
  the chip name from the annotation file header as the platform name.

  \code{load.affy.mapping} creates and populates a new mapping, using
  the version of the NetAffx release as the mapping name.
}
\value{
  For \code{read.affy.anno}, a data frame representing a subset of
  annotations suitable for importing into GT.DB, with key=value
  attributes from the file header stored as attributes of the result.

  For \code{load.affy.platform} and \code{load.affy.mapping}, the result
  is the number of assays processed.
}
\seealso{
  \code{\link{load.affy.chp.data}},
  \code{\link{mk.platform}}, \code{\link{mk.assay}},
  \code{\link{mk.mapping}}, \code{\link{mk.assay.position}}.
}
\examples{
\dontrun{
anno <- read.affy.anno('Axiom_GW_Hu_SNP.r2.na30.annot.csv')
load.affy.platform(anno, 'Affymetrix Axiom 600K')
load.affy.mapping(anno)
}
}
\keyword{IO}
