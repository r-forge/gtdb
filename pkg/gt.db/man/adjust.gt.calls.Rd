%
% Copyright (C) 2009, Perlegen Inc.
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
\name{adjust.gt.calls}
\alias{adjust.gt.calls}
\title{Manually Assign Genotype Calls}
\description{
  Interactively assign genotype calls by inspection of a plot of
  reference versus alternate signal intensities.
}
\usage{
adjust.gt.calls(data, ..., radius=6)
}
\arguments{
  \item{data}{a data frame returned by \code{\link{reshape.gt.data}}.}
  \item{\dots}{additional arguments passed to
    \code{\link{gt.cluster.plot}}.}
  \item{radius}{the radius around the mouse pointer selected by
    clicking, in points.}
}
\details{
  The intensity data is first plotted, then the user is iteratively
  prompted for a genotype value (in a/h/b/n notation) in the console
  window.  Once a value has been entered, the plot window becomes
  active and the mouse pointer can be used to select points whose
  genotypes should be changed to that value.  Each left click causes a
  zone around the mouse pointer position to be updated.  The process
  is terminated by right clicking and selecting \sQuote{Stop} from the
  menu.  The user is then prompted for another genotype value in the
  console window.  To finish here, enter \code{q}.
}
\value{
  A data frame like \code{data} but with updated genotype calls.  If a
  column called \code{orig.genotype} does not already exist, then it
  will be created and populated with the original unmodified genotypes.
}
\seealso{
  \code{\link{fetch.gt.data}},
  \code{\link{reshape.gt.data}},
  \code{\link{gt.cluster.plot}}.
}
\examples{
\dontrun{
gt <- fetch.gt.data('Demo_2',raw.data=TRUE)
d <- reshape.gt.data(gt[32,], na.codes='n')
adjust.gt.calls(d, scales=list(log=TRUE))
}
}
\keyword{manip}
\keyword{iplot}
