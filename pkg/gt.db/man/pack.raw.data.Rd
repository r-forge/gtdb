%
% Copyright (C) 2010, 23andMe, Inc.
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
\name{pack.raw.data}
\alias{pack.raw.data}
\alias{unpack.raw.data}
\title{Pack and Unpack Raw Genotype Data}
\description{
  Convert raw genotype data between multi-column data frames and packed
  raw vectors.
}
\usage{
pack.raw.data(data, raw.layout)
unpack.raw.data(raw, raw.layout)
}
\arguments{
  \item{data}{a data frame.}
  \item{raw}{a raw vector.}
  \item{raw.layout}{one of \code{'signal'}, \code{'seqread'}, or
    \code{'chpdata'}.}
}
\details{
  In a packed raw vector, all the data elements for an individual sample
  are grouped together, to simplify subsetting by sample.  For
  \code{pack.raw.data}, columns will be coerced to the appropriate data
  types before packing.

  The raw data layouts have the following structures:

  \describe{
    
    \item{\code{'signal'}}{two columns: \code{signal.a} and
      \code{signal.b}.  These represent transformed intensities for the
      A and B alleles.  These should fit into 16-bit unsigned integers
      (0 to 65535, with 65535 reserved for representing \code{NA}).}
    
    \item{\code{'seqread'}}{four columns: \code{fwd.a}, \code{rev.a},
      \code{fwd.b}, and \code{rev.b}.  These represent forward and
      reverse sequencing read counts for the A and B alleles.  Each
      should fit into a single unsigned byte, with the value of 255
      reserved for \code{NA}.}
    
    \item{\code{'chpdata'}}{three columns: \code{log.ratio},
      \code{strength}, and \code{forced.call}.  These represent the
      base-2 log of the intensity ratio for A and B alleles; the base-2
      log of the overall signal strength; and a forced genotype call
      with no quality thresholding.  The first two are stored as
      fixed-point signed 8.8-bit numbers in the raw vector form.}

  }
}
\value{
  For \code{pack.raw.data}, a raw vector encoding the source data
  according to \code{raw.layout}.

  For \code{unpack.raw.data}, a data frame representing the contents of
  the raw vector interpreted according to \code{raw.layout}.
}
\seealso{
  \code{\link{mk.assay.data}}, \code{\link{reshape.gt.data}}.
}
\examples{
x <- data.frame(signal.a=c(100,10), signal.b=c(10,100))
y <- pack.raw.data(x, 'signal')
unpack.raw.data(y, 'signal')

x <- data.frame(fwd.a=c(2,5), rev.a=c(4,4), fwd.b=c(3,0), rev.b=c(4,0))
y <- pack.raw.data(x, 'seqread')
unpack.raw.data(y, 'seqread')

x <- data.frame(log.ratio=c(0.1,3.5), strength=c(5.5,7.5),
                forced.call=c('h','a'))
y <- pack.raw.data(x, 'chpdata')
unpack.raw.data(y, 'chpdata')
}
\keyword{manip}
