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
\name{use.gt.db}
\alias{use.gt.db}
\title{Define GT.DB Database Connection}
\description{
  Define a DBI connection to be used to access a GT.DB database.  This
  connection is used implicitly by other GT.DB functions that interact
  with the database.
}
\usage{
use.gt.db(dbConnection)
}
\arguments{
  \item{dbConnection}{a connection object from
    \code{\link[DBI:dbConnect]{dbConnect}}.}
}
\value{
  None.
}
\seealso{
  \code{\link[DBI:dbDriver]{dbDriver}},
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link{init.gt.db}}.
}
\examples{\dontrun{
library(RSQLite)
tmpname <- tempfile('db')
gt.db <- dbConnect(dbDriver('SQLite'), tmpname)
use.gt.db(gt.db)
}}
\keyword{database}
