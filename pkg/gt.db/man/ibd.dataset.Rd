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
  Pairwise identity-by-descent (IBD) is estimated using data from a
  series of intervals of size \code{binsz}, equally spaced across the
  genome.  The algorithm determines the minimum number of alleles
  identical by state (IBS) for SNPs within each bin, and averages
  these values across bins.  The minimum IBS value gives an upper
  limit on IBD for that bin, and approaches IBD as the number of
  assayed SNPs increases.

  Determining IBD for close relatives requires only a fraction of the
  available data from a whole genome scan.  Very accurate estimates of
  IBD are not particularly useful, because of intrinsic variability in
  the recombination process.  The default values for \code{binsz},
  \code{part}, and \code{parts} should not need to be changed.

  Bins with few SNPs may not be sufficiently informative for IBD, if
  there is a substantial probability of unrelated samples sharing the
  same haplotype by chance.  Rare genotyping errors can also impact
  apparent IBD, since a single error within a bin will change the
  minimum IBS.  The \code{min.snps} setting ensures reasonable
  informativeness, and \code{max.snps}, \code{maf.min},
  \code{hw.p.min}, and \code{gt.rate.min} control the expected error
  rate.

  The estimated genomic proportion with IBD=1 is biased upwards
  because a small proportion of bins are not sufficiently informative
  to distinguish IBD=0 from IBD=1.  The estimate for IBD=2 seems to
  be less biased.  The following table shows expected genomic IBD
  proportions for common familial relationships.

  \tabular{cccl}{
   IBD=0 \tab IBD=1 \tab IBD=2 \tab Relationship \cr
   1.00 \tab 0.00 \tab 0.00 \tab unrelated \cr
   0.75 \tab 0.25 \tab 0.00 \tab first cousins \cr
   0.50 \tab 0.50 \tab 0.00 \tab half siblings \cr
   0.00 \tab 1.00 \tab 0.00 \tab parent-child \cr
   0.25 \tab 0.50 \tab 0.25 \tab full siblings \cr
   0.00 \tab 0.00 \tab 1.00 \tab duplicate samples \cr
  }

  It is important to note that the actual IBD proportions for most
  relative pairs (except for parent-child) show substantial natural
  variation, independent of the variation due to the estimation
  procedure.
}
\value{
  A list with three elements.  The first two elements are matrices,
  containing the estimated proportions of the genome with IBD=1 and
  IBD=2, for all pairs of samples in the dataset.  The third element
  is an ordered list of sample names from the GTPUB dataset.
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
