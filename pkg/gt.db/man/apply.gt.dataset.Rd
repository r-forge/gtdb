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
\name{apply.gt.dataset}
\alias{apply.gt.dataset}
\title{Apply a Function to a Genotype Dataset}
\description{
  Apply a function to a genotype dataset, processed in manageable 
  chunks, and aggregate the results.
}
\usage{
apply.gt.dataset(gt.dataset, part.fn, aggr.fn, ..., GT.DATA='gt.data')
}
\arguments{
  \item{gt.dataset}{a dataset description from \code{\link{gt.dataset}}.}
  \item{part.fn}{a function to apply to each part of the dataset.}
  \item{aggr.fn}{an aggregator function to accumulate results.}
  \item{\dots}{additional arguments passed to \code{part.fn}.}
  \item{GT.DATA}{the name to use for passing the genotype argument to
     \code{part.fn}.}
}
\details{
  Based on the specification in \code{gt.dataset}, the dataset is
  loaded in chunks defined by \code{parts} \code{binsz}, and
  \code{gt.filter}.  \code{part.fn} is applied iteratively to each
  chunk.  The genotype data is passed as the \code{GT.DATA} argument
  to \code{part.fn}.  The aggregator function \code{aggr.fn} should
  accept two arguments: an accumulated result, initially set to NULL;
  and a result for a single chunk.
}
\value{
  Depends on form of \code{part.fn} and \code{aggr.fn}.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{gt.dataset}}.
}
\examples{
gt.demo.check()
s <- apply.gt.dataset(gt.dataset('Demo_1'), summary,
                      rbind, GT.DATA='object')
str(s)
s <- apply.gt.dataset(gt.dataset('Demo_1'), summary,
                      aggr.fn=function(a,b) c(a,list(b)),
                      GT.DATA='object')
str(s)
}
\keyword{database}
