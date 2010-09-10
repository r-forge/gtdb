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
\name{store.prcomp}
\alias{store.prcomp}
\alias{rm.prcomp}
\title{Store or Remove Principal Components Results}
\description{
  Store, or remove, principal components analysis results in the database.
}
\usage{
store.prcomp(x, prcomp.name, description,
             nc=ncol(x$loadings), is.hidden=FALSE)
rm.prcomp(dataset.name, prcomp.name)
}
\arguments{
  \item{x}{a principal components result returned by
    \code{\link{prcomp.gt.data}} or \code{\link{prcomp.gt.dataset}}.}
  \item{prcomp.name}{a short unique identifier for the analysis.}
  \item{description}{a free-text description of this analysis.}
  \item{nc}{the number of components to store.}
  \item{is.hidden}{logical: indicates if the dataset is hidden.}
  \item{dataset.name}{the name of a parent dataset for which an analysis
    is to be removed.}
}
\details{
  Analysis names only need to be unique within the scope of a
  particular dataset.
}
\value{
  The number of principal components stored into the database.
}
\seealso{
  \code{\link{fetch.prcomp}},
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
