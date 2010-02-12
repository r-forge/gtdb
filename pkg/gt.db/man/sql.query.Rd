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
\name{sql.query}
\alias{sql.query}
\alias{sql.exec}
\alias{sql.exec-methods}
\alias{sql.exec,MySQLConnection-method}
\alias{sql.exec,OraConnection-method}
\alias{sql.exec,SQLiteConnection-method}
\title{Simplified SQL Statement Execution}
\description{
  These functions parse and execute SQL statements, and either fetch
  query results or return the number of affected rows.  All result
  sets are also cleaned up.  Some database differences and differences
  in DBI implementations are also concealed.
}
\usage{
  sql.query(db, sql, ...)
  sql.exec(db, sql, ..., chunk.kb=256, progress=FALSE)
}
\arguments{
  \item{db}{a DBI connection object returned by
    \code{\link[DBI:dbConnect]{dbConnect}}.}
  \item{sql}{a valid SQL statement, optionally with embedded
    bind variables.}
  \item{\dots}{optionally, a data frame or arguments to be used to
    construct a data frame of bind variables.}
  \item{chunk.kb}{when processing multiple rows of bind variables,
    a rough limit on the amount of data to send per query.}
  \item{progress}{logical: indicates whether to report progress during
    long operations.}
}
\details{
  To facilitate database agnostic code, several special elements can
  be used in SQL statements, which are internally replaced with
  database specific forms:

  \describe{
    \item{:user:}{the database user name.}
    \item{:sysdate:}{the current date.}
    \item{:unhex:(...)}{a function to convert from a hex string to binary.}
    \item{:clob:(...)}{a wrapper for passing long character data back to R.}
    \item{:blob:(...)}{a wrapper for passing long binary data back to
      R, rendered as a hexadecimal string.}
    \item{:fromdual:}{a placeholder for the Oracle \code{FROM DUAL}
      idiom for \code{SELECT} statements that do not need table data.}
  }

  Both \code{sql.query} and \code{sql.exec} emulate prepared
  statements and bind variables in a consistent way across databases.
  The emulation uses the Oracle syntax for bind variables (i.e.,
  \sQuote{:1}, \sQuote{:2}, etc in SQL are substituted with values of
  the bind variables).

  A single call to these functions may result in multiple SQL
  statements.  These will be grouped into a single transaction for
  databases that support that, so that a call to \code{sql.exec}
  should be \dQuote{all or none}.  RMySQL currently does not support
  transactions, so a call to \code{sql.exec} may partially fail.
}
\value{
  For \code{sql.query}, a data frame constructed from the results
  returned by the query.  Column names are transformed by converting
  \samp{\_} to \samp{.} and changing to lower case.  For \code{sql.exec},
  the return value is the number of affected rows.
}
\seealso{
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link[DBI:fetch]{fetch}}.
}
\examples{
\dontrun{
sql.exec(db, 'create table xyzzy (a number, b number)')
sql.exec(db, 'insert into xyzzy values (:1,:2)', a=6, b=1:4)
sql.exec(db, 'insert into xyzzy values (:1,:2)', data.frame(1,2))
sql.query(db, 'select * from xyzzy')
sql.exec(db, 'drop table xyzzy')
sql.exec(db, 'select :user:, :sysdate: :fromdual:')
}
}
\keyword{database}

