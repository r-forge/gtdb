%
% Copyright (C) 2009, Perlegen Sciences, Inc.
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
\name{init.gt.db}
\alias{init.gt.db}
\title{Initialize GT.DB Database}
\description{
  Create all the standard GT.DB database objects (tables, indexes) in an
  empty database.
}
\usage{
init.gt.db(db.mode=c('raw','hex','zip'))
}
\arguments{
  \item{db.mode}{storage mode for packed objects in the database: either
    \code{'hex'}, \code{'raw'}, or \code{'zip'}.  See details.}
}
\details{
  This should be called after connecting to a new database using
  \code{\link[DBI:dbConnect]{dbConnect}} and \code{\link{use.gt.db}}.
  Scripts for creating GT.DB tables and indexes are installed under
  \code{library(help='gt.db')$path} in the \file{schema} subdirectory.

  The \code{db.mode} argument controls how genotypes, quality scores,
  and raw data are stored in the database.  The default (\code{'raw'})
  is to store genotypes as character data, and quality scores and
  underlying data as binary blobs.  In \code{'hex'} mode, quality scores
  and underlying data are stored as strings of hex digits.  In
  \code{'zip'} mode, all are stored in compressed form.

  Without appropriate plugins, the SQLite interface supports only
  \code{'hex'} mode.  At the moment, \code{'zip'} mode is only supported
  in MySQL.
}
\seealso{
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link{use.gt.db}},
  \code{\link{demo}(setup.gt.demo)}.
}
\examples{\dontrun{
# create temporary in-memory SQLite database
dbx <- dbConnect(dbDriver('SQLite'), ':memory:')
use.gt.db(dbx)
init.gt.db(db.mode='hex')
demo('setup.gt.demo')
}}
\keyword{database}
