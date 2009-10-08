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
\name{ls.test}
\alias{ls.test}
\title{List Association Test Result Sets}
\description{
  This returns a list of association test result sets in the database
  for a specified genotype dataset.
}
\usage{
ls.test(dataset.name, test.name='\%', show.all=FALSE, show.ids=FALSE)
}
\arguments{
  \item{dataset.name}{a genotype dataset name.}
  \item{test.name}{an SQL \code{LIKE} expression for matching
    test set names.}
  \item{show.all}{logical: indicates if hidden result sets should be
    included in the output.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per test set, and 8 or 9 columns:
  \item{test.id}{if \code{show.ids} is set: the unique integer key
    for this test set.}
  \item{test.name}{a short identifier for the test set, unique when
    combined with \code{term} for the specified genotype dataset.}
  \item{description}{a free-text description of the test set.}
  \item{fit}{the type of model fit: \code{'lm'}, \code{'glm'}, etc.}
  \item{model}{the model formula used for scoring each genotype assay.}
  \item{term}{used for differentiating among multipe test result sets
    associated with a single analysis}
  \item{is.hidden}{logical: indicates if the test set is hidden.}
  \item{created.by}{the user name that created the test set.}
  \item{created.dt}{the creation date of the test set.}
}
\seealso{
  \code{\link{fetch.test.scores}}, \code{\link{store.test.scores}},
  \code{\link{rm.test}}.
}
\examples{
# FIXME
}
\keyword{database}
