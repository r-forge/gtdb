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
\name{fetch.prcomp}
\alias{fetch.prcomp}
\title{Load Principal Components Results for a Genotyping Dataset}
\description{
  Load principal components analysis results for a Genotyping Dataset.
}
\usage{
fetch.prcomp(dataset.name, prcomp.name, nc)
}
\arguments{
  \item{dataset.name}{the unique identifier for the dataset.}
  \item{prcomp.name}{the unique identifier for the analysis.}
  \item{nc}{the number of components to retrieve.}
}
\details{
  If \code{prcomp.name} is not present, then if there is just one
  visible analysis for this dataset, that one will be loaded.
}
\value{
  The same sort of structure as returned by \code{\link{prcomp.gt.data}}.
}
\seealso{
  \code{\link{store.prcomp}},
  \code{\link{prcomp.gt.data}},
  \code{\link{prcomp.gt.dataset}}.
}
\examples{
\dontrun{
pc <- prcomp(gt.dataset('Demo_1', gt.filter=(ploidy=='A')))
store.prcomp(pc, 'demo_pc_1', 'Demo PCA results')
fetch.prcomp('Demo', 'demo_pc_1', nc=4)
}
}
\keyword{database}
