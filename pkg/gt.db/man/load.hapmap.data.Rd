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
\name{load.hapmap.data}
\alias{load.hapmap.data}
\title{Import HapMap Genotype Data}
\description{
  Import non-redundant, forward-stranded genotype data from the 
  International HapMap Project.
}
\usage{
load.hapmap.data(files, project.name='HapMap', map=TRUE, verbose=TRUE)
}
\arguments{
  \item{files}{a vector of phased or unphased HapMap genotype data
    file names.}
  \item{project.name}{the project to be associated with this data.}
  \item{map}{logical: indicates whether to load assay and map
    information into the database, or to assume it is already there.}
  \item{verbose}{logical: indicates whether to report progress.}
}
\details{
  This function supports loading recent HapMap Phase II and Phase III
  genotype data.  If files include data from several population panels,
  then the genotype data is merged by rsID into a single dataset
  spanning all those panels.  Only non-redundant forward-orientation
  files are supported.

  It has been tested against Phase II unphased r22 and r24, Phase II
  phased r22, and Phase III unphased r2 and r3.

  In principle, it should be possible to define a \dQuote{Phase II}
  platform once, and have multiple HapMap datasets refer to that
  platform, with release-specific map information.  In practice, this is
  challenging because releases can mix and match data from multiple
  assays of the same SNP (i.e. CEU on one assay, JPT+CHB and YRI on
  another).  We do not currently load the underlying assay information
  and instead treat each rsID in a release as an \dQuote{assay}.

  It is possible to load both phased and unphased datasets in the same
  release so they refer to the same mapping.  First load the unphased
  dataset, then load the phased data with \code{map=FALSE}.
}
\seealso{
  \code{\link{hapmap.subjects}}.
}
\examples{
\dontrun{
# create HapMap project
demo('setup.hapmap')

# get file list for latest Phase II HapMap release
base <- 'http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes'
release  <- '/latest_phaseII_ncbi_b36/fwd_strand/non-redundant'
path <- paste(base, release, sep='')
re <- '.*(genotypes_chr[^_]+_..._[0-9A-Za-z_.]+gz).*'
files <- sub(re, '\\1', grep(re, readLines(path), value=TRUE))

# download to current directory
for (f in files) {
    download.file(paste(path,f,sep='/'), f)
}

# load into current database
load.hapmap.data(files)
}
}
\keyword{IO}
