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
\name{mk.sample}
\alias{mk.sample}
\alias{rm.sample}
\title{Create or Remove Dataset Samples}
\description{
  \code{mk.sample} and \code{rm.sample} insert or remove samples
  from a genotype dataset.
}
\usage{
mk.sample(dataset.name, data)
rm.sample(dataset.name, sample.name)
}
\arguments{
  \item{dataset.name}{a short unique identifier for the dataset.}
  \item{data}{a data frame with one row per sample.  See details.}
  \item{sample.name}{a vector of sample identifiers.}
}
\details{
  A data frame of sample information should have (at least) four
  columns: \code{sample.name}, \code{subject.name}, \code{gender}, and
  \code{position}.  This information is important because it affects how
  various functions interpret genotyping data in this dataset.
  Sample gender is used to determine ploidy for sex linked assays, and
  the position is used to index arrays of packed genotype data.
  Positions are 1-based integers.
}
\value{
If successful, the number of rows inserted or deleted from the sample
table.
}
\seealso{
  \code{\link{ls.sample}}, \code{\link{mk.subject}}.
}
\examples{\dontrun{
s <- data.frame(
    sample.name='NA_12345_1'
    subject.name='NA_12345',
    gender='M',
    position=1
)
mk.sample('Demo_1', s)
rm.sample('Demo_1, s$sample.name)
}}
\keyword{database}
