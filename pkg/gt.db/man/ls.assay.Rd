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
\name{ls.assay}
\alias{ls.assay}
\title{List Assay Definitions}
\description{
  List definitions of assays associated with a genotyping platform.
}
\usage{
ls.assay(platform.name, show.ids=FALSE)
}
\arguments{
  \item{platform.name}{a genotyping platform name.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per assay, and up to 6 columns:
  \item{assay.id}{if \code{show.ids} is set: the unique integer ID for
    this assay.}
  \item{assay.name}{the assay name.}
  \item{flags}{an integer value composed of single-bit flags.}
  \item{alleles}{a slash separated list of allele sequences.}
  \item{probe.seq}{the genomic flanking sequence for the assay, with
    the variant position denoted by an underscore (\code{'_'}).}
  \item{alt.name}{an alternate name for the assay, if available.}
}
\seealso{
  \code{\link{ls.assay.position}}, \code{\link{mk.assay}}.
}
\examples{
gt.demo.check()
head(ls.assay('Demo_Set_1'))
head(ls.mapping('Demo_Set_1', show.ids=TRUE))
}
\keyword{database}
