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
\name{ls.dataset}
\alias{ls.dataset}
\title{List Genotype Datasets}
\description{
  This returns a list of genotype datasets defined in the current database.
}
\usage{
ls.dataset(project.name='\%', dataset.name='\%',
           show.all=FALSE, show.ids=FALSE)
}
\arguments{
  \item{project.name}{an SQL \code{LIKE} expression for matching
    project names.}
  \item{dataset.name}{an SQL \code{LIKE} expression for matching
    dataset names.}
  \item{show.all}{logical: indicates if hidden datasets and members of
    hidden projects should be included in the output.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per dataset, and 8 or 11 columns:
  \item{dataset.id}{if \code{show.ids} is set: the unique integer key
    for this dataset.}
  \item{dataset.name}{a short, unique identifier for the dataset.}
  \item{project.id}{if \code{show.ids} is set: the unique integer key
    for the project associated with this dataset.}
  \item{project.name}{the project associated with this dataset.}
  \item{platform.id}{if \code{show.ids} is set: the unique integer key
    for the platform associated with this dataset.}
  \item{platform.name}{the genotyping platform associated with this
    dataset.}
  \item{description}{a free-text description of the dataset.}
  \item{raw.layout}{a keyword indicating how raw data associated with
    individual assays is structured.}
  \item{is.hidden}{logical: indicates if the project is hidden.}
  \item{created.by}{the user name that created the dataset.}
  \item{created.dt}{the creation date of the dataset.}
}
\seealso{
  \code{\link{ls.project}}, \code{\link{ls.platform}},
  \code{\link{mk.dataset}}.
}
\examples{
gt.demo.check()
ls.dataset()
ls.dataset(show.ids=TRUE)
}
\keyword{database}
