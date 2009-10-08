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
\name{fetch.test.scores}
\alias{fetch.test.scores}
\title{Load Association Test Results for a Genotyping Dataset}
\description{
  Load association test results for a Genotyping Dataset.
}
\usage{
fetch.test.scores(dataset.name, test.name, term, max.pvalue)
}
\arguments{
  \item{dataset.name}{the unique identifier for a dataset.}
  \item{test.name}{the unique identifier for an analysis.}
  \item{term}{if present, for result sets with multiple results per
    assay, specifies which set of results to fetch.}
  \item{max.pvalue}{if present, specifies an upper bound for P values to
    be fetched.}
}
\value{
  The same sort of structure as returned by \code{\link{score.gt.data}}.
}
\seealso{
  \code{\link{store.test.scores}},
  \code{\link{score.gt.data}}
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
