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
\name{gt.split}
\alias{gt.split}
\alias{gt.paste}
\title{Convert between Packed Genotype Strings and Genotype Vectors}
\description{
  \code{gt.split} takes a string of a/h/b/n genotypes and returns a
  vector with one genotype per element; \code{gt.paste} reverses
  this conversion.
}
\usage{
gt.split(s, convert=c('score.b','score.a','char','none'),
         alleles, na.codes=c(), strand=c('+','-'), sep='/')
gt.paste(v, convert=c('score.b','score.a','char','none'),
         alleles, na.codes=c(), strand=c('+','-'), sep='/')
}
\arguments{
  \item{s}{a character string of packed a/h/b/n genotypes.}
  \item{v}{a vector of unpacked genotypes.  See details.}
  \item{convert}{how genotypes should be represented.}
  \item{alleles}{a vector of length 2 specifying A and B alleles,
    for \code{convert='char'}.}
  \item{na.codes}{a vector of missing genotype codes: see details.}.
  \item{strand}{for \code{convert='char'}, strand to report.}
  \item{sep}{for \code{convert='char'}, a string for separating
    alleles.}
}
\details{
  If \code{convert='score.a'}, then unpacked genotypes are coded as
  numeric scores, where \code{"a"}=2, \code{"h"}=1, and \code{"b"}=0.
  If \code{convert='score.b'}, then unpacked genotypes are coded as
  numeric scores, where \code{"a"}=0, \code{"h"}=1, and \code{"b"}=2.
  if \code{convert='char'}, then genotypes are coded as factors with
  levels formed from concatenated pairs of alleles.  If
  \code{convert='none'}, then genotypes are coded as factors with
  three levels, \code{'a'}, \code{'h'}, and \code{'b'}.

  Missing genotypes are by default coded as \code{NA} in the unpacked
  format.  If \code{na.codes} is not empty, then the specified single
  character codes are instead reported as separate factor levels in
  the unpacked format.
}
\value{
  For \code{gt.split}, a vector of genotypes coded as specified by
  \code{convert}.  For \code{gt.paste}, a packed genotype string.
}
\seealso{
  \code{\link{fetch.gt.data}}, \code{\link{unpack.gt.matrix}}.
}
\examples{
gt <- 'aahabbnnaxa'
gt.split(gt)
gt.split(gt, convert='char', alleles=c('A','C'))
gt.split(gt, convert='char', alleles=c('A','C'), strand='-')
gt.split(gt, convert='char', alleles=c('A','C'),
         na.codes='n', sep='')
gt.split(gt, convert='none', na.codes='n')
gt.split(gt, convert='none', na.codes=c('n','x'))

x <- gt.split(gt, convert='char', alleles=c('A','C'))
gt.paste(x, convert='char', alleles=c('A','C'))
gt.paste(x, convert='char', alleles=c('C','A'))
gt.paste(c(0,1,NA,2,1,2,NA), convert='score.b')
}
\keyword{manip}
