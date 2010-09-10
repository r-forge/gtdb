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
\name{store.subject.data}
\alias{store.subject.data}
\title{Store Subject Data}
\description{
  Stores subject phenotype data into the database from a data frame.
}
\usage{
store.subject.data(project.name, data)
}
\arguments{
  \item{project.name}{the short unique identifier for the project.}
  \item{data}{a data frame of subject phenotypes.  See details.}
}
\details{
  The \code{data} argument should have a \code{subject.name} column,
  with subjects previously defined by \code{\link{mk.subject}}.
  Additional columns should supply values of subject attributes
  previously defined by \code{\link{mk.subject.attr}}.  Missing values
  are not explicitly recorded in the database.
}
\value{
  A vector of counts of rows inserted in the database for each column
  of phenotype information in the input data.
}
\seealso{
  \code{\link{mk.subject}},
  \code{\link{mk.subject.attr}},
  \code{\link{fetch.subject.data}},
  \code{\link{fetch.pt.data}}.
}
\examples{
\dontrun{
p <- ls.subject('Demo')['subject.name']
p$stuff <- rnorm(nrow(p))
mk.subject.attr('Demo', p['stuff'], 'Some random stuff')
store.subject.data('Demo', p)
head(fetch.subject.data('Demo')
rm.subject.attr('Demo', 'stuff')
}
}
\keyword{database}
