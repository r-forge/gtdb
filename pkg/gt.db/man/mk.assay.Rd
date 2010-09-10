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
\name{mk.assay}
\alias{mk.assay}
\title{Create Genotyping Assay Definitions}
\description{
  \code{mk.assay} defines genotyping assays associated with an already
  defined genotyping platform.
}
\usage{
mk.assay(platform.name, data, progress=FALSE)
}
\arguments{
  \item{platform.name}{a short unique identifier for the platform.}
  \item{data}{a data frame with one row per assay.  See details.}
  \item{progress}{logical: indicates whether to report progress while
    the data is loaded.}
}
\details{
  A data frame of assay information can provide up to five columns:
  \describe{
    \item{assay.name}{a unique name (within this platform) for the assay.}
    \item{flags}{an optional integer value composed of single-bit flags.}
    \item{alleles}{a slash separated list of valid alleles.}
    \item{probe.seq}{an optional genomic flanking sequence for the assay,
      with the variant position denoted by an underscore (\code{'_'}).}
    \item{alt.name}{an alternate name for the assay, which must be
      unique within this platform, if present.}
  }
}
\value{
  If successful, the number of rows inserted into the assay table.
}
\seealso{
  \code{\link{ls.assay}}, \code{\link{mk.assay.position}}.
}
\examples{\dontrun{
data(demo_01)
head(assay.def.01)
mk.assay('Demo_Set_1', assay.def.01)
}}
\keyword{database}
