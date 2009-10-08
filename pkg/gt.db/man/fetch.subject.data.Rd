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
\name{fetch.subject.data}
\alias{fetch.subject.data}
\title{Load Subject Data for a Genotyping Project}
\description{
  Returns a data frame of subject phenotypes for the specified
  project, mapped to the appropriate R datatypes.
}
\usage{
fetch.subject.data(project.name, cols, show.all=FALSE)
}
\arguments{
  \item{project.name}{the short unique identifier for the project.}
  \item{cols}{a vector of column names to be included.}
  \item{show.all}{logical: indicates if hidden attributes should be
    included in the output.}
}
\value{
  A data frame with one row per subject, and one column per attribute
  defined for this project.  The column datatypes are determined using
  the information in the \code{SUBJECT_ATTR} table.  The row order is
  arbitrary, but row names are set to the sample names.
}
\seealso{
  \code{\link{mk.subject}},
  \code{\link{ls.subject}},
  \code{\link{store.subject.data}}.
}
\examples{
gt.demo.check()
head(fetch.subject.data('Demo'))
}
\keyword{database}
