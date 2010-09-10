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
\name{ibd.plot}
\alias{ibd.plot}
\title{Plot Identity by Descent Data}
\description{
  Plot identity by descent data for all sample pairs.
}
\usage{
ibd.plot(ibd, jitter.amount, split=0.3,
         panel.width=list(0.8,'npc'))
}
\arguments{
  \item{ibd}{an IBD estimate from \code{\link{ibd.dataset}}.}
  \item{jitter.amount}{scale of random noise to be added to the IBD=2
    values to resolve overlapping observations.}
  \item{split}{the proportion of plot height devoted to the IBD=1
    histogram.}
  \item{panel.width}{the proportion of plot width to use for plotting
    data, as opposed to axes, tick labels, etc.}
}
\details{
  Produces a two-panel plot, with the upper panel showing pairwise
  IBD=1 versus IBD=2 proportions, and the lower panel showing a
  histogram of the IBD=1 values.  In the upper panel, outliers are
  colored by whether they are outliers on IBD=1, on IBD=2, or both.
  Outliers on IBD=1 are defined by a Bonferroni corrected P value of
  less than 0.05, assuming a normal distribution for unrelated
  individuals.
}
\value{
  A plot object of class \code{"trellis"}.
}
\seealso{
  \code{\link{ibd.dataset}},
  \code{\link{ibd.summary}},
  \code{\link{ibd.outliers}}.
}
\examples{
gt.demo.check()
ibd <- ibd.dataset('Demo_1', min.snps=10, binsz=10e6, parts=1)
ibd.plot(ibd, jitter=0.005)
}
\keyword{hplot}
