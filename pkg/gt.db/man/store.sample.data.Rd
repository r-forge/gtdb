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
\name{store.sample.data}
\alias{store.sample.data}
\title{Store Sample Data}
\description{
  Stores sample phenotype data into the database from a data frame.
}
\usage{
store.sample.data(dataset.name, data)
}
\arguments{
  \item{dataset.name}{the short unique identifier for the dataset.}
  \item{data}{a data frame of sample phenotypes.  See details.}
}
\details{
  The \code{data} argument should have a \code{sample.name} column,
  with samples previously defined by \code{\link{mk.sample}}.
  Additional columns should supply values of sample attributes
  previously defined by \code{\link{mk.sample.attr}}.  Missing values
  are not explicitly recorded in the database.
}
\value{
  A vector of counts of rows inserted in the database for each column
  of phenotype information in the input data.
}
\seealso{
  \code{\link{mk.sample}},
  \code{\link{mk.sample.attr}},
  \code{\link{fetch.sample.data}},
  \code{\link{fetch.pt.data}}.
}
\examples{
\dontrun{
p <- ls.sample('Demo_1')['sample.name']
p$stuff <- rnorm(nrow(p))
mk.sample.attr('Demo_1', p['stuff'], 'Some random stuff')
store.sample.data('Demo_1', p)
head(fetch.sample.data('Demo_1')
rm.sample.attr('Demo_1', 'stuff')
}
}
\keyword{database}
