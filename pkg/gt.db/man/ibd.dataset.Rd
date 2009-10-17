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
\name{ibd.dataset}
\alias{ibd.dataset}
\title{Calculate Identity by Descent for a Genotype Dataset}
\description{
  Estimates identity by descent for all sample pairs in a genotype
  dataset, from observed identity by state information.
}
\usage{
ibd.dataset(dataset.name, binsz=1e6, part=1, parts=10,
            maf.min=0.1, hw.p.min=0.01, gt.rate.min=0.98, ...)
}
\arguments{
  \item{dataset.name}{the identifier for the dataset to analyze.}
  \item{binsz, part, parts}{passed to \code{\link{fetch.gt.data}},
    to identify the regions to be analyzed.}
  \item{maf.min}{minimum minor allele frequency to be included in the 
    analysis.}
  \item{hw.p.min}{minimum Hardy-Weinberg equilibrium P value for SNPs
    to be included in the analysis.}
  \item{gt.rate.min}{minimum call rate for SNPs to be included.}
  \item{\dots}{additional arguments for \code{\link{ibd.gt.data}}.}
}
\details{
  This is a wrapper around \code{\link{ibd.gt.data}} that helps to
  select a suitable genome-wide autosomal subset of a genotype dataset
  for IBD analysis.  The \code{maf.min}, \code{hw.p.min}, and
  \code{gt.rate.min} filters improve the IBD estimates by excluding
  problematic data that may have elevated error rates.
}
\value{
  A list of two square matrices with row and column names set to sample
  names in the input data, containing the estimated proportions of the
  genome with IBD=1 and IBD=2, for all pairs of samples in the dataset.
}
\seealso{
  \code{\link{ibd.plot}},
  \code{\link{ibd.summary}},
  \code{\link{ibd.outliers}}.
}
\examples{
gt.demo.check()
ibd <- ibd.dataset('Demo_1', min.snps=10, binsz=10e6, parts=1)
ibd.plot(ibd, jitter=0.005)
}
\keyword{manip}
