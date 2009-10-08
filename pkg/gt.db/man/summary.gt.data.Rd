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
\name{summary.gt.data}
\alias{summary.gt.data}
\alias{summary.gt.dataset}
\title{Genotype Data Summary}
\description{
  Computes genotype counts, allele frequencies, call rates, and tests
  for Hardy Weinberg equilibrium for a packed genotype data structure.
}
\usage{
\method{summary}{gt.data}(object, sample.mask, by.sample=FALSE, ...)
\method{summary}{gt.dataset}(object, sample.mask, by.sample=FALSE, ...)
}
\arguments{
  \item{object}{either a data frame of genotypes from
    \code{\link{fetch.gt.data}}, or a dataset description from
    \code{\link{gt.dataset}}.} 
  \item{sample.mask}{an optional mask identifying a subset of samples
    to be included in the summaries.}
  \item{by.sample}{logical: indicates if statistics should be computed
    by sample (as opposed to by SNP).}
  \item{\dots}{not used.}
}
\details{
  The \code{sample.mask} argument may either be a character string or
  a logical vector with one element per genotype.
}
\value{
  If \code{by.sample} is \code{FALSE}: a data frame with one row per
  genotype assay, and 9 columns:

  \item{NN}{count of missing (\code{'n'}) genotypes.}
  \item{AA, AB, BB}{counts of diploid r/h/a genotypes.}
  \item{A_, B_}{counts of haploid r/a genotypes.}
  \item{gt.rate}{genotype call rate for this SNP.}
  \item{freq.a, freq.b}{allele frequencies.}
  \item{hw.p.value}{Hardy-Weinberg equilibrium P value, calculated
                    from a likelihood ratio test.}

  If \code{by.sample} is \code{TRUE}: a data frame with one row per
  sample and 9 columns, with samples ordered by position in the
  genotype strings, and row names set to \code{RNA_DNA_SOURCE_ID}.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{gt.dataset}}.
}
\examples{
gt.demo.check()
g <- fetch.gt.data('Demo_1')[1:10,]
summary(g)
summary(g, attr(g,'gender')[['F']])
summary(g, by.sample=TRUE)[1:10,]
}
\keyword{manip}
