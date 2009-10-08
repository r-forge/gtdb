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
\name{mk.assay.group}
\alias{mk.assay.group}
\alias{rm.assay.group}
\title{Create or Remove a Genotype Assay Group}
\description{
  These functions insert or remove entries from the genotyping
  assay group table. 
}
\usage{
mk.assay.group(platform.name, assay.group, description)
rm.assay.group(platform.name, assay.group)
}
\arguments{
  \item{platform.name}{the unique platform identifier.}
  \item{assay.group}{an identifier for the assay group.}
  \item{description}{a free-text description of the assay group.}
}
\details{
  An assay group and its constituent assays can only be removed if it
  is not in use.  Any mappings or datasets that refer to these assays
  must be explicitly deleted first.
}
\value{
  If successful, the number of rows inserted in or removed from the
  assay group table (i.e., typically, 1). 
}
\seealso{
  \code{\link{ls.assay.group}},
  \code{\link{ls.platform}},
  \code{\link{mk.platform}},
  \code{\link{mk.assay}}.
}
\examples{\dontrun{
mk.platform('My_600K', 'My 600K Array Set')
mk.assay.group('My_600K', 'Des01', '300K Array 1 of 2')
mk.assay.group('My_600K', 'Des02', '300K Array 2 of 2')
rm.assay.group('My_600K', c('Des01','Des02'))
rm.platform('My_600K')
}}
\keyword{database}
