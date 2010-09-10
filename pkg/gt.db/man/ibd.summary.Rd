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
\name{ibd.summary}
\alias{ibd.summary}
\alias{ibd.outliers}
\title{Summarize Identity by Descent Analysis Results}
\description{
  Summarize results from an identity-by-descent analysis.
}
\usage{
ibd.summary(ibd, alpha=0.05)
ibd.outliers(ibd, alpha=0.05)
}
\arguments{
  \item{ibd}{results from \code{\link{ibd.dataset}}.}
  \item{alpha}{a P value threshold for identifying IBD=1 outliers.}
}
\details{
  \code{ibd.summary} and \code{ibd.outliers} both identify sample
  pairs with elevated genome-wide IBD=1 or IBD=2 proportions.  For
  IBD=1, a threshold is selected assuming a normal distribution for
  unrelated individuals, using the specified \code{alpha} threshold
  with a Bonferroni correction.  For IBD=2, an (arbitrary) threshold
  of 0.05 is used.

  Unrelated pairs should have neither elevated; parent-sibling pairs,
  cousins and half sibs should have only IBD=1 elevated; duplicates
  should have only IBD=2 elevated; siblings should have both elevated.
}
\value{
  For \code{\link{ibd.summary}}, a table of counts of distinct sample
  pairs, indicating how many have neither IBD=1 or IBD=2 elevated;
  only elevated IBD=1; only elevated IBD=2; or both elevated.

  For \code{\link{ibd.outliers}}, a data frame with 7 columns:
  \item{rds.id.1, sample.name.1}{identifiers for sample 1.}
  \item{rds.id.2, sample.name.2}{identifiers for sample 2.}
  \item{grp}{the type of outlier: one of \code{'IBD=1'},
    \code{'IBD=2'}, \code{'Both'}.}
  \item{ibd.1}{the inferred proportion of the genome with IBD=1.}
  \item{ibd.2}{the inferred proportion of the genome with IBD=2.}
}
\seealso{
  \code{\link{ibd.dataset}},
  \code{\link{ibd.plot}}.
}
\examples{
gt.demo.check()
ibd <- ibd.dataset('Demo_1', min.snps=10, binsz=10e6, parts=1)
ibd.summary(ibd)
ibd.outliers(ibd)
}
\keyword{manip}
