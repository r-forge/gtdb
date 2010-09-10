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
\name{fetch.gt.data}
\alias{fetch.gt.data}
\title{Load Genotype Data for a Genotyping Dataset}
\description{
  Returns a data frame of assay information and genotypes for the
  specified dataset.
}
\usage{
fetch.gt.data(dataset.name, mapping.name, assay.name, dbsnp.rsid,
              part, parts, by=c('assay','position'), binsz=1, where,
              show.ids=FALSE, genotype=TRUE, qscore=FALSE, raw.data=FALSE)
}
\arguments{
  \item{dataset.name}{the unique identifier for the dataset.}
  \item{mapping.name}{an identifier for the assay mapping to use.}
  \item{assay.name}{a vector of assay names.}
  \item{dbsnp.rsid}{a vector of integer dbSNP rsID values.}
  \item{part}{a vector of subsets 1..parts.}
  \item{parts}{the number of subsets to split the dataset into.}
  \item{by}{specifies how to construct subsets.  See details.}
  \item{binsz}{specifies how to construct subsets.  See details.}
  \item{where}{additional SQL \code{WHERE} clauses to limit the data
    returned.}
  \item{show.ids}{logical: indicates whether to include values of
    database keys.}
  \item{genotype}{logical: indicates if genotypes should be fetched.}
  \item{qscore}{logical: indicates if quality scores should be fetched.}
  \item{raw.data}{logical: indicates if underlying raw data should be
    fetched, if available.}
}
\details{
  If an assay mapping is not explicitly specified, and there is just
  one visible mapping for this dataset's platform, that one will be
  used by default.

  The \code{assay.name}, \code{dbsnp.rsid}, and \code{part} arguments
  are mutually exclusive ways of specifying a set of SNPs for which
  data should be retrieved.

  If \code{part} and \code{parts} are specified, then the dataset is
  divided into \code{parts} roughly similarly sized subsets, and
  \code{part} specifies which subset to retrieve.  The chunking can be
  done either by primary key (\code{by='assay.id'}) or by position
  (\code{by='position'}).  To preserve some locality, assay subsets
  are interleaved such that \code{trunc(by/binsz)\%\%parts==(part-1)}.
}
\value{
  A data frame with one row per SNP.  If any flags have been defined
  for this dataset, each flag is expanded to a logical column.
  \item{assay.data.id}{the unique integer identifier for this entry in
    the \code{ASSAY_DATA} table.}
  \item{assay.id}{if \code{show.ids} is set: the unique key for this
    assay.}
  \item{assay.name}{the name of this assay.}
  \item{alleles}{a slash-separated list of alleles for this assay.}
  \item{scaffold}{an identifier for the scaffold to which this assay
    has been mapped.}
  \item{position}{the 1-based genomic position in the specified scaffold.}
  \item{strand}{the strand (+/-) to which the assay was mapped.}
  \item{ploidy}{one of \code{"A"} (autosomal), \code{"X"} (X linked),
    \code{"Y"} (Y linked), or \code{"M"} (mitochondrial).}
  \item{dbsnp.rsid}{the dbSNP rsID to which this assay is mapped.}
  \item{dbsnp.orient}{the relative orientation of the dbSNP refsnp
    cluster.}
  \item{genotype}{a string of a/h/b/n/x genotypes.}
  \item{qscore}{a string of packed quality scores.}
  \item{raw.data}{a string containing packed raw data (such as
    signal intensities or read counts).}
}
\seealso{
  \code{\link{summary.gt.data}}, \code{\link{unpack.gt.matrix}}.
}
\examples{
gt.demo.check()
gt <- fetch.gt.data('Demo_1')
str(gt)
nrow(fetch.gt.data('Demo_1', part=1, parts=4))
nrow(fetch.gt.data('Demo_1', part=1, parts=100))
summary.gt.data(gt[1:10,])
}
\keyword{database}
