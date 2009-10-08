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
\name{store.test.scores}
\alias{store.test.scores}
\alias{rm.test}
\title{Store or Remove Association Test Results}
\description{
  Store, or remove, association test results in the database.
}
\usage{
store.test.scores(x, test.name, description, is.hidden=FALSE)
rm.test(dataset.name, test.name)
}
\arguments{
  \item{x}{a principal components result returned by
    \code{\link{prcomp.gt.data}} or \code{\link{prcomp.gt.dataset}}.}
  \item{test.name}{a short unique identifier for the analysis.}
  \item{description}{a free-text description of this analysis.}
  \item{is.hidden}{logical: indicates if the result set is hidden.}
  \item{dataset.name}{the name of a parent dataset for which an analysis
    is to be removed.}
}
\details{
  Analysis names only need to be unique within the scope of a
  particular dataset.
}
\value{
  The number of results stored into the database.
}
\seealso{
  \code{\link{fetch.test.scores}},
  \code{\link{score.gt.data}}.
}
\examples{
\dontrun{
pt <- fetch.pt.data('Demo_2')
gt <- fetch.gt.data('Demo_2')
pt$status <- as.logical(rbinom(nrow(pt), 1, 0.5))
x <- score.gt.data(status~genotype, pt, gt)
store.test.scores(x, 'status_1', 'Test analysis')
y <- fetch.test.scores('Demo_2', 'status_1')
str(y)
}
}
\keyword{database}
