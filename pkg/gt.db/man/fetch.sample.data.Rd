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
\name{fetch.sample.data}
\alias{fetch.sample.data}
\title{Load Sample Data for a Genotype Dataset}
\description{
  Returns a data frame of sample phenotypes for the specified dataset,
  mapped to the appropriate R datatypes.
}
\usage{
fetch.sample.data(dataset.name, cols, show.all=FALSE)
}
\arguments{>
  \item{dataset.name}{the short unique identifier for the dataset.}
  \item{cols}{a vector of column names to be included.}
  \item{show.all}{logical: indicates if hidden attributes should be
    included in the output.}
}
\value{
  A data frame with one row per sample, and one column per attribute
  defined for this dataset.  The column datatypes are determined using
  the information in the \code{SAMPLE_ATTR} table.  The row order is
  arbitrary, but row names are set to the sample names.  The result
  also includes the predefined sample attributes: subject name,
  gender, and dataset position.
}
\seealso{
  \code{\link{mk.sample}}, \code{\link{ls.sample}},
  \code{\link{fetch.pt.data}}, \code{\link{store.sample.data}}.
}
\examples{
gt.demo.check()
head(fetch.sample.data('Demo_1'))
}
\keyword{database}
