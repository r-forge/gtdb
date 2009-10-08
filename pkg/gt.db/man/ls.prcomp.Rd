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
\name{ls.prcomp}
\alias{ls.prcomp}
\title{List Principal Components Result Sets}
\description{
  This returns a list of principal components analysis result sets in
  the database for a specified genotype dataset.
}
\usage{
ls.prcomp(dataset.name, prcomp.name='\%', show.all=FALSE, show.ids=FALSE)
}
\arguments{
  \item{dataset.name}{a genotype dataset name.}
  \item{prcomp.name}{an SQL \code{LIKE} expression for matching
    principal components analysis set names.}
  \item{show.all}{logical: indicates if hidden result sets should be
    included in the output.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per result set, and 9 or 10 columns:
  \item{prcomp.id}{if \code{show.ids} is set: the unique integer key
    for this test set.}
  \item{prcomp.name}{a short identifier for the result set, which is
    unique for the specified genotype dataset.}
  \item{description}{a free-text description of the result set.}
  \item{fn.call}{the text of the call used to perform the analysis.}
  \item{components}{the number of components retained in the database.}
  \item{samples}{the number of samples included in the analysis.}
  \item{assays}{the number of assays included in the analysis.}
  \item{is.hidden}{logical: indicates if the result set is hidden.}
  \item{created.by}{the user name that created the result set.}
  \item{created.dt}{the creation date of the result set.}
}
\seealso{
  \code{\link{fetch.prcomp}}, \code{\link{store.prcomp}},
  \code{\link{rm.prcomp}}.
}
\examples{
# FIXME
}
\keyword{database}
