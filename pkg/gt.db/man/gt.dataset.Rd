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
\name{gt.dataset}
\alias{gt.dataset}
\title{Genotype Dataset Specification}
\description{
  Describes a genotype dataset in the database and how to process it.
}
\usage{
gt.dataset(dataset.name, gt.filter=TRUE, parts=10,
           by='assay', binsz=1, progress=interactive())
}
\arguments{
  \item{dataset.name}{the unique identifier for the dataset.}
  \item{gt.filter}{an expression to use for subsetting the genotype
    data.}
  \item{parts, by, binsz}{specify how to construct chunks.
    See \code{\link{fetch.gt.data}}.}
  \item{progress}{logical: indicates whether to display a progress
    bar when processing the dataset.}
}
\details{
  The specified dataset is loaded in chunks defined by \code{parts}
  and \code{binsz}, and \code{part.fn} is applied iteratively to each
  chunk.  The genotype data is passed as the \code{gt.data} argument
  to \code{part.fn}.  The aggregator function \code{aggr.fn} should
  accept two arguments: an accumulated result, initially set to NULL;
  and a result for a single chunk.
}
\value{
  A list of class \code{gt.dataset}, with elements representing all
  the arguments of \code{gt.dataset}.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{apply.gt.dataset}}.
}
\examples{
gt.demo.check()
s <- apply.gt.dataset(gt.dataset('Demo_1'), summary,
                      rbind, GT.DATA='object')
str(s)
s <- apply.gt.dataset(gt.dataset('Demo_1', gt.filter=(ploidy=='A')),
                      summary, rbind, GT.DATA='object')
str(s)
}
\keyword{database}
