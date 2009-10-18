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
\name{init.gt.db}
\alias{init.gt.db}
\title{Initialize GT.DB Database}
\description{
  Create all the standard GT.DB database objects (tables, indexes) in an
  empty database.
}
\usage{
init.gt.db()
}
\details{
  This should be called after connecting to a new database using
  \code{\link[DBI:dbConnect]{dbConnect}} and \code{\link{use.gt.db}}.
  Scripts for creating GT.DB tables and indexes are installed under
  \code{library(help='gt.db')$path} in the \file{schema} subdirectory.
}
\seealso{
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link{use.gt.db}},
  \code{\link{demo}(setup.gt.demo)}.
}
\examples{\dontrun{
# create new SQLite database in a temporary file
fn <- tempfile()
dbx <- dbConnect(dbDriver('SQLite'), fn)
# unlink it so it will go away at the end of the session
unlink(fn)
use.gt.db(dbx)
init.gt.db()
demo('setup.gt.demo')
}}
\keyword{database}
