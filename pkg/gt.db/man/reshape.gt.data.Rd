%
% Copyright (C) 2009, Perlegen Sciences, Inc.
% Copyright (C) 2010, 23andMe, Inc.
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
\name{reshape.gt.data}
\alias{reshape.gt.data}
\title{Reshape Genotype Data}
\description{
  This takes a data frame of genotypes with one row per assay, and
  returns a \dQuote{long} data frame with one row per genotype per
  sample.
}
\usage{
reshape.gt.data(gt.data, ...)
}
\arguments{
  \item{gt.data}{a data frame of genotypes from \code{fetch.gt.data}.}
  \item{\dots}{additional arguments passed to \code{\link{gt.split}}.}
}
\value{
  A data frame with one row per assay per sample, and columns
  corresponding to available data extracted from \code{gt.data}.  This
  normally includes genotypes, and may include quality scores and/or
  one or more columns of additional underlying raw data.
  \item{assay.name}{an assay identifier.}
  \item{sample.name}{a sample identifier.}
  \item{genotype}{the genotype for the specified sample and assay.}
  \item{qscore}{the corresponding quality score.}
  \item{\dots}{additional columns of raw data, depending on the
    dataset.  For datasets with signal intensities, there will be two
    columns: \code{signal.a} and \code{signal.b}.  For data with
    sequencing read counts, there will be four columns: \code{fwd.a},
    \code{rev.a}, \code{fwd.b}, and \code{rev.b}.  For Affymetrix CHP
    data, there will be three columns: \code{log.ratio},
    \code{strength}, and \code{forced.call}.}
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{unpack.gt.matrix}},
  \code{\link{gt.split}}, \code{\link{unpack.raw.data}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_2',raw.data=TRUE)
head(reshape.gt.data(gt))
d <- reshape.gt.data(gt[seq(32,232,40),], na.codes='n')
gt.cluster.plot(d)
}
\keyword{manip}
