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
\name{load.affy.chp.data}
\alias{load.affy.chp.data}
\title{Import Affymetrix CHP Genotype Data}
\description{
  Import Affymetrix CHP genotype data files into the database, including
  genotypes, quality scores, and underlying measurements.
}
\usage{
load.affy.chp.data(dataset, anno, files, progress=TRUE)
}
\arguments{
  \item{dataset}{the unique identifier of the dataset to receive the
    genotype data.}
  \item{anno}{an annotation data frame from \code{\link{read.affy.anno}}.}
  \item{files}{a character vector of all the text CHP file names to be
    imported.}
  \item{progress}{logical: indicates whether to report progress during
    the database load.}
}
\details{
  The input files should be tab-delimited text versions of CHP files as
  created by the \code{apt-chp-to-txt} utility from the Affymetrix Power
  Tools.

  Loading CHP data is a two stage process.  In the first stage, we
  reorganize the input data into new files that represent all samples
  across smaller sets of SNPs.  In the second stage, we load data from
  these files into the database.  The two stages effectively allow us to
  perform an \dQuote{out-of-core} transpose of the genotype matrix,
  because we typically will not be able to hold an entire genotype
  dataset for a large study in memory.  The temporary files require
  about 20\% of the space of the original text CHP files.  They are
  created in the current directory and deleted at the end of the import.

  When creating a dataset to receive CHP data, be sure to specify
  \code{raw.format='chpdata'}.
}
\value{
  The number of assays for which data was loaded.
}
\seealso{
  \code{\link{read.affy.anno}},
  \code{\link{mk.dataset}},
  \code{\link{mk.assay.data}}.
}
\keyword{IO}
