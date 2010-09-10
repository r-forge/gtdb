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
\name{fetch.pt.data}
\alias{fetch.pt.data}
\title{Load Phenotype Data for a Genotype Dataset}
\description{
  Returns a data frame of merged sample and subject phenotypes for the
  specified genotype dataset, mapped to the appropriate R datatypes.
}
\usage{
fetch.pt.data(dataset.name, cols, pca=FALSE, show.all=FALSE)
}
\arguments{
  \item{dataset.name}{the short unique identifier for the dataset.}
  \item{cols}{a vector of either sample or subject attributes to be
    included.} 
  \item{pca}{a logical indicating whether the current principal
    components analysis results should be included, or a list of
    arguments to be passed to \code{\link{fetch.prcomp}} to specify
    what to load.}
  \item{show.all}{logical: indicates if hidden attributes should be
    included in the output.}
}
\value{
  A data frame with one row per sample, and one column per attribute
  defined for this dataset, combining sample data from
  \code{\link{fetch.sample.data}} and subject data from
  \code{\link{fetch.subject.data}}.  The column datatypes are
  determined using the information in the \code{SAMPLE_ATTR} and
  \code{SUBJECT_ATTR} tables.  The row order conforms to the order of
  samples in the genotype data, and row names are set to the sample
  names.  The result also includes the predefined sample attributes:
  subject name, gender, and dataset position.
}
\seealso{
  \code{\link{fetch.sample.data}}.
  \code{\link{fetch.subject.data}}.
  \code{\link{fetch.prcomp}}.
}
\examples{
gt.demo.check()
head(fetch.pt.data('Demo_1'))
head(fetch.pt.data('Demo_1', cols='plate'))
}
\keyword{database}
