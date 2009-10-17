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
\name{ibd.gt.data}
\alias{ibd.gt.data}
\title{Estimate Pairwise Identity by Descent for Genotype Data}
\description{
  Estimate identity-by-descent for all sample pairs in a set of
  genotype data.
}
\usage{
ibd.gt.data(gt.data, binsz=1e6, min.snps=25, max.snps=50,
            min.gt=0.8, ibs.limit=2)
}
\arguments{
  \item{gt.data}{a data frame of genotypes from \code{\link{fetch.gt.data}}.}
  \item{binsz}{Size (in base pairs) of bins to use for estimating IBD.}
  \item{min.snps}{The minimum number of SNPs to consider when estimating
    IBD.  Bins with fewer SNPs will be excluded.}
  \item{max.snps}{The maximim number of SNPs to consider.  Additional
    SNPs in a bin will be ignored.}
  \item{min.gt}{minimum sample call rate for inclusion in the analysis.}
  \item{ibs.limit}{bins with average IBS \code{ibs.limit} fold higher
    than the median IBS are excluded due to low information content.}
}
\details{
  Pairwise identity-by-descent (IBD) is estimated by subdividing the
  input data into intervals of intervals of size \code{binsz}, equally
  spaced across the genome.  The algorithm determines the minimum number
  of alleles identical by state (IBS) for SNPs within each bin, and
  averages these values across bins.  The minimum IBS value gives an
  upper limit on IBD for that bin, and approaches IBD as the number of
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
  informativeness, and \code{max.snps} helps to limit the error rate.

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
gt <- fetch.gt.data('Demo_2')
ibd <- ibd.gt.data(subset(gt, ploidy=='A'))
}
\keyword{manip}
