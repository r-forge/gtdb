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
\name{ls.assay.position}
\alias{ls.assay.position}
\title{List Assay Positions}
\description{
  List mapped genomic positions of assays associated with a genotyping
  platform.
}
\usage{
ls.assay.position(platform.name, mapping.name, show.ids=FALSE)
}
\arguments{
  \item{platform.name}{a genotyping platform name.}
  \item{mapping.name}{an identifier for the assay mapping to use.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\details{
  If \code{mapping.name} is missing, it will default to the current
  (visible) mapping for the specified platform, if that is unique.
}
\value{
  A data frame with one row per assay, and 7 or 8 columns:
  \item{assay.id}{if \code{show.ids} is set: the unique integer ID for
    this assay.}
  \item{assay.name}{the assay name.}
  \item{scaffold}{a string identifying the sequence to which the
    assay is mapped.}
  \item{position}{a one-based position within the specified scaffold.}
  \item{strand}{either \code{"+"}, \code{"-"}, or \code{NA}.}
  \item{ploidy}{describes the expected allele count: see
    \code{\link{mk.assay.position}}.}
  \item{dbsnp.rsid}{the dbSNP refSNP cluster ID for this assay.}
  \item{dbsnp.orient}{the orientation of the assay compared to the
    dbSNP cluster: either \code{"+"}, \code{"-"}, or \code{NA}.}
}
\seealso{
  \code{\link{ls.assay}}, \code{\link{mk.assay.position}}.
}
\examples{
gt.demo.check()
head(ls.assay.position('Demo_Set_1'))
}
\keyword{database}
