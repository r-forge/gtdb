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
\name{ls.sample}
\alias{ls.sample}
\title{List Samples in a Genotype Dataset}
\description{
  This returns a description of all samples defined in the specified
  genotype dataset.
}
\usage{
ls.sample(dataset.name, show.ids=FALSE)
}
\arguments{
  \item{dataset.name}{the unique identifier for the dataset.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
}
\value{
  A data frame with one row per dataset, and 4 or 6 columns, sorted by
  \code{position}:
  \item{sample.id}{if \code{show.ids} is set: the unique integer key
    for this sample.}
  \item{sample.name}{a short, unique identifier for the sample.}
  \item{subject.id}{if \code{show.ids} is set: the unique integer key
    for the subject associated with this sample.}
  \item{subject.name}{a short, unique identifier for the subject.}
  \item{gender}{a factor indicating the gender of this sample.}
  \item{position}{a 1-based integer indicating the position of this
    sample's genotyping data in the genotype arrays for this dataset.}
}
\seealso{
  \code{\link{ls.dataset}}, \code{\link{mk.sample}}.
}
\examples{
gt.demo.check()
head(ls.sample('Demo_1'))
head(ls.sample('Demo_1', show.ids=TRUE))
}
\keyword{database}
