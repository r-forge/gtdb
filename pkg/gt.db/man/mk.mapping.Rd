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
\name{mk.mapping}
\alias{mk.mapping}
\alias{rm.mapping}
\title{Create or Remove a Mapping for a Genotyping Platform}
\description{
  These functions insert or remove entries from the platform mapping
  table.
}
\usage{
mk.mapping(platform.name, mapping.name, description,
           assembly, is.hidden=FALSE)
rm.mapping(platform.name, mapping.name)
}
\arguments{
  \item{platform.name}{the unique platform identifier.}
  \item{mapping.name}{an identifier for the mapping.}
  \item{description}{a free-text description of the mapping.}
  \item{assembly}{the target assembly for the mapping.  This is used to
    determine whether two mappings have compatible coordinates.}
  \item{is.hidden}{logical: indicates if the mapping is hidden.}
}
\details{
  A platform may have multiple mapping sets, for different assemblies,
  or even for the same assembly with different alignment parameters.
  Removing a mapping results in deletion of all the associated assay
  map positions.
}
\value{
  If successful, the number of rows inserted in or removed from the
  mapping table (i.e., typically, 1). 
}
\seealso{
  \code{\link{ls.mapping}},
  \code{\link{mk.platform}},
  \code{\link{mk.assay.position}}.
}
\examples{\dontrun{
mk.platform('My_600K', 'My 600K Array Set')
mk.mapping('My_600K', 'b36_1', 'Build 36 Vender Annotations', 'ncbi_b36')
rm.mapping('My_600K', 'b36_1')
rm.platform('My_600K')
}}
\keyword{database}
