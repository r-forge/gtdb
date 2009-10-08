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
\name{nw.align}
\alias{nw.align}
\alias{nw.score}
\title{Needleman and Wunsch Sequence Alignment}
\description{
  Returns a global or partial-local optimal alignment of two nucleotide
  sequences, with linear gap penalty.  \code{nw.score} returns just
  the alignment score.
}
\usage{
nw.align(gt1, gt2, gap=-1, pm=1, mm=0, ends=FALSE)
nw.score(gt1, gt2, gap=-1, pm=1, mm=0, ends=FALSE)
}
\arguments{
  \item{gt1, gt2}{nucleotide sequences to be aligned.}
  \item{gap}{gap penalty.}
  \item{pm}{score for a perfect match.}
  \item{mm}{score for a mismatch.}
  \item{ends}{logical: specifies whether overhanging ends are not
    penalized.  See details.}
}
\details{
  The default is to perform a global alignment across the full lengths
  of both sequences.  The \code{ends} argument can be used to specify
  that overhangs at one or both ends of one or both sequences should
  not be penalized.  If \code{ends} is a scalar, then it controls all
  four possible overhangs.  It may also be set to a vector of four
  logical values, to control scoring of each possible overhang: the
  left and right sides of \code{gt1}, and left and right sides of
  \code{gt2}.

  The dynamic programming algorithm in \code{nw.align} is O(m*n) in
  time and space, and hence is not appropriate for aligning very long
  sequences.  The method in \code{nw.score} is O(n) in space.

  To get the longest common subsequence, use 
}
\value{
  A list:
  \item{score}{the optimal alignment score.}
  \item{pct.matched}{over the aligned intervals of the two sequences,
    the percentage identity.}
  \item{ends}{the starting and ending base positions of the aligned
    intervals in \code{gt1} and \code{gt2}.}
  \item{alignment}{a vector of three strings: the aligned \code{gt1},
    a representation of the matched positions, and the aligned \code{gt2}.}
}
\references{
  Needleman, S. B., & Wunsch, C.D. (1970) A general method applicable
  to the search for similarities in the amino acid sequence of two
  proteins.  \emph{J. Mol. Biol.} \bold{48}: 443-453.
}
\seealso{
  \code{\link{nw.orient.assay}}.
}
\examples{
s1 <-       'GAGTTC'
s2 <- 'AACATCGAGGTCTACGGATT'
nw.align(s1, s2, ends=FALSE)
nw.align(s1, s2, ends=TRUE)
}
\keyword{manip}
