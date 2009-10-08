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
\name{score.and.store}
\alias{score.and.store}
\title{Test SNPs for Association and Store Results}
\description{
  Perform a series of single-point SNP association tests on a genotype
  dataset, using an arbitrary scoring function, and store the results
  in the \code{TEST} and \code{TEST_RESULT} tables.
}
\usage{
score.and.store(dataset.name, test.name, description, formula,
                score.fn=NULL, pt.filter=TRUE, gt.filter=TRUE,
                pca=FALSE, part=1:parts, parts=200,
                dryrun=FALSE, ...)
}
\arguments{
  \item{dataset.name}{the name of the dataset to be analyzed.}
  \item{test.name}{a unique identifier to associate with these
    test results.} 
  \item{description}{a description of this test result set.}
  \item{formula}{a symbolic description of the model to be fit.}
  \item{score.fn}{the scoring function to be applied to each SNP.}
  \item{pt.filter}{an expression to use for subsetting on the
    phenotype table.}
  \item{gt.filter}{an expression to use for subsetting on the
    genotype table.}
  \item{pca}{a logical indicating whether the current principal
    components analysis results for this dataset should be loaded, or
    a list of arguments to be passed to \code{\link{fetch.prcomp}} to
    specify an analysis to be loaded.}
  \item{part}{a vector of slices 1..parts.}
  \item{parts}{the number of slices to split the dataset into.}
  \item{dryrun}{logical: if \code{TRUE}, then do not write results to
    the database.}
  \item{\dots}{additional arguments to pass to \code{\link{score.gt.data}}.}
}
\details{
  This is essentially a wrapper around \code{\link{score.gt.data}},
  that handles fetching phenotype and genotype data, and storing
  results back into the database.  The \code{formula},
  \code{score.fn}, \code{pt.filter}, and \code{gt.filter} arguments
  are passed to \code{\link{score.gt.data}}.  The \code{pca} argument
  is passed to \code{\link{fetch.pt.data}}.  And the \code{part} and
  \code{parts} arguments are used with \code{\link{fetch.gt.data}}.
  The \code{description} argument is used for constructing new entries
  in the \code{TEST} table in GTPUB.

  If executed as part of an LSF \dQuote{job array}, this function will
  automatically initialize \code{parts} on each node to divide the
  computation evenly across the job array.
}
\seealso{
  \code{\link{fetch.pt.data}},
  \code{\link{fetch.gt.data}},
  \code{\link{score.gt.data}}.
}
\examples{
gt.demo.check()
## score 5\% of the SNPs in Demo_1 dataset
score.and.store('Demo_1', 'Test_1', 'Example Analysis',
                plate ~ genotype, score.chisq,
                gt.filter=(pass==1),
                part=1, parts=20, dryrun=TRUE,
                progress=TRUE)
}
\keyword{htest}
\keyword{regression}
\keyword{database}
