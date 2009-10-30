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
\name{mk.dataset}
\alias{mk.dataset}
\alias{rm.dataset}
\title{Create or Remove Genotype Datasets}
\description{
  These functions insert or delete entries from the genotype dataset
  table.
}
\usage{
mk.dataset(dataset.name, project.name, platform.name,
           description, raw.layout=c(NA,'signal','seqread'),
           is.hidden=FALSE)
rm.dataset(dataset.name)
}
\arguments{
  \item{dataset.name}{a short unique identifier for the dataset.}
  \item{project.name}{the project associated with this dataset.}
  \item{platform.name}{the genotyping platform for this dataset.}
  \item{description}{a free-text description of the dataset.}
  \item{raw.layout}{a keyword describing how raw data associated with
    individual assays is structured.}
  \item{is.hidden}{logical: indicates if the dataset is hidden.}
}
\details{
  While datasets are nested within projects, dataset names must be
  globally unique.
}
\value{
If successful, the number of rows inserted in or deleted from the
dataset table (i.e., 1). 
}
\seealso{
  \code{\link{ls.dataset}}, \code{\link{set.hidden}}.
}
\examples{
\dontrun{
mk.dataset('Demo_9', project.name='Demo', platform.name='Demo',
           description='Demo Dataset #9')
rm.dataset('Demo_9')
}
}
\keyword{database}
