#
# Copyright (C) 2009, Perlegen Sciences, Inc.
# Copyright (C) 2010, 23andMe, Inc.
#
# Written by David A. Hinds <dhinds@sonic.net>
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the license, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#---------------------------------------------------------------------

# Functions for hiding some SQL messiness and database differences

#---------------------------------------------------------------------

# prepare arguments to be bound to a prepared statement

.sql.prep.data <- function (...)
{
    data <- as.data.frame(lapply(list(...), unname))
    if (ncol(data) > 0) {
        fix <- (sapply(data, class) == "factor")
        data[, fix] <- lapply(data[, fix, drop=FALSE], as.character)
        fix <- (sapply(data, class) == "logical")
        data[, fix] <- lapply(data[, fix, drop=FALSE], as.numeric)
    }
    data
}

# clear all pending result sets

.dbClearAll <- function(db)
{
    invisible(lapply(dbListResults(db), dbClearResult))
}

#---------------------------------------------------------------------

# The following code is used to emulate "prepared statements" that
# work the same way as they do for Oracle.  While this doesn't yield
# the same performance advantage, it does simplify query writing.
#
# Use caution: some of these functions take shortcuts for parsing
# SQL statements with regular expressions and may mangle valid SQL.

.myPrepareStatement <- function(db, sql, ...)
{
    # Caution: this regular expression ignores quoting
    sql <- gsub('%', '%%', sql)
    list(db=db, sql=gsub(':([0-9+])\\b', "%\\1$s", sql))
}

.myBindStatement <- function(stmt, data)
{
    quote <- function(x)
    {
        e <- gsub("'", "''", x)
        if.na(x, 'NULL', paste("'", e, "'", sep=''))
    }
    do.call(sprintf, c(stmt$sql, lapply(data, quote)))
}

.myExecStatement <- function(stmt, data, ...)
{
    dbSendQuery(stmt$db, .myBindStatement(stmt, data[1,,drop=FALSE]), ...)
}

#---------------------------------------------------------------------

# the following functions are used to hide some SQL dialect issues,
# and work around the fact that ROracle does not support CLOB/BLOBs

.split.clob <- function(max.len=16000)
{
    start <- seq(1, max.len, 4000)
    paste('dbms_lob.substr(\\1,4000,', start, ') \\1_',
          seq_along(start), collapse=', ', sep='')
}
.split.blob <- function(fn.1='', fn.2='', max.len=32000)
{
    start <- seq(1, max.len, 2000)
    paste(fn.1, '(dbms_lob.substr(', fn.2, '(\\1),2000,', start,
          ')) \\1_', seq_along(start), collapse=', ', sep='')
}

.fixup.sql <- function(db, sql)
{
    if (class(db) == 'OraConnection') {
        sql <- gsub(':user:', 'user', sql)
        sql <- gsub(':sysdate:', 'sysdate', sql)
        sql <- gsub(':clob:\\((\\w+)\\)', .split.clob(), sql)
        sql <- gsub(':blob:\\((\\w+)\\)', .split.blob(), sql)
        sql <- gsub(':hex.blob:\\((\\w+)\\)',
                    .split.blob('rawtohex'), sql)
        sql <- gsub(':unhex:', 'hextoraw', sql)
        sql <- gsub(':fromdual:', 'from dual', sql)
        sql <- gsub(':floor:', 'floor', sql)
        sql <- gsub(':mod:', 'mod', sql)
        sql <- gsub('%%', ',', sql)
    }
    if (class(db) == 'MySQLConnection') {
        sql <- gsub(':user:', 'user()', sql)
        sql <- gsub(':sysdate:', 'sysdate()', sql)
        sql <- gsub(':clob:\\((\\w+)\\)', '\\1', sql)
        sql <- gsub(':unzip.clob:\\((\\w+)\\)', 'uncompress(\\1) \\1', sql)
        sql <- gsub(':blob:\\((\\w+)\\)', '\\1', sql)
        sql <- gsub(':hex.blob:\\((\\w+)\\)', 'hex(\\1) \\1', sql)
        sql <- gsub(':unhex:', 'unhex', sql)
        sql <- gsub(':zip:', 'compress', sql)
        sql <- gsub(':zip.unhex:\\(([^)]+)\\)', 'compress(unhex(\\1))', sql)
        sql <- gsub(':hex.unzip.blob:\\((\\w+)\\)',
                    'hex(uncompress(\\1)) \\1', sql)
        sql <- gsub(':fromdual:', '', sql)
        sql <- gsub(':floor:', 'floor', sql)
        sql <- gsub(':mod:', 'mod', sql)
        sql <- gsub('%%', ',', sql)
    }
    if (class(db) == 'SQLiteConnection') {
        sql <- gsub(':user:', 'null', sql)
        sql <- gsub(':sysdate:', 'current_timestamp', sql)
        sql <- gsub(':clob:\\((\\w+)\\)', '\\1', sql)
        sql <- gsub(':blob:\\((\\w+)\\)', '\\1', sql)
        sql <- gsub(':hex.blob:\\((\\w+)\\)', 'hex(\\1) \\1', sql)
        sql <- gsub(':unhex:', 'unhex', sql)
        sql <- gsub(':fromdual:', '', sql)
        sql <- gsub(':floor:\\(', 'round(-0.5+', sql)
        sql <- gsub(':mod:', '', sql)
        sql <- gsub('%%', ')%(', sql)
    }
    sql
}

.fixup.result  <- function(db, sql, result)
{
    pack <- function(...) paste(..., sep='')
    grps <- function(x)
    {
        m <- gregexpr(pack(x,'\\(\\w+\\)'), sql)[[1]]
        if (m[1] < 0) return(c())
        mapply(substr, sql, m+7, m-2+attr(m,'match.length'))
    }
    if (class(db) == 'OraConnection') {
        for (n in c(grps(':blob:'), grps(':clob:'))) {
            cols <- grep(sprintf('^%s_[0-9]+',n), names(result),
                         ignore.case=TRUE)
            result[,n] <- do.call(pack, lapply(result[cols],if.na,''))
            result <- result[-cols]
        }
    }
    names(result) <- chartr('A-Z_', 'a-z.', names(result))
    result
}

#---------------------------------------------------------------------

sql.query <- function(db, sql, ...)
{
    efn <- function(e) {
        .dbClearAll(db)
        stop(e)
    }
    data <- .sql.prep.data(...)
    sql.2 <- .fixup.sql(db, sql)

    if (class(db) == 'OraConnection') {
        ps <- tryCatch(dbPrepareStatement(db, sql.2, data), error=efn)
        rs <- tryCatch(dbExecStatement(ps, data), error = efn)
    } else {
        ps <- tryCatch(.myPrepareStatement(db, sql.2, data), error=efn)
        rs <- tryCatch(.myExecStatement(ps, data), error = efn)
    }

    r <- tryCatch(fetch(rs, n = -1), error = efn)
    dbClearResult(rs)
    .fixup.result(db, sql, r)
}

#---------------------------------------------------------------------

.chunk.rows <- function(data, chunk.kb)
{
    rowsz <- object.size(data) / nrow(data)
    max(1, floor(1024 * chunk.kb / rowsz))
}

.ora.sql.exec <- function(db, sql, ..., chunk.kb=256, progress=FALSE)
{
    efn <- function(e)
    {
        dbRollback(db)
        .dbClearAll(db)
        stop(e)
    }
    data <- .sql.prep.data(...)
    sql <- .fixup.sql(db, sql)

    chunk.rows <- .chunk.rows(data, chunk.kb)
    ps <- tryCatch(dbPrepareStatement(db, sql, data), error=efn)
    nr <- 0
    nd <- nrow(data)
    if (nd && progress) progress.bar(0, nd)
    for (lo in seq(1,max(1,nd),chunk.rows)) {
        d <- data[lo:min(lo+chunk.rows-1,nd),,drop=FALSE]
        tryCatch(dbExecStatement(ps, d), error=efn)
        nr <- nr + nrow(d)
        if (nd && progress) progress.bar(nr, nd)
    }

    nr <- dbGetRowsAffected(ps)
    dbCommit(db)
    dbClearResult(ps)
    nr
}

# MySQL can insert multiple rows with a single SQL statement, and this
# is much more efficient than doing multiple independent inserts.

.myInsertMultiple <- function(stmt, data, ...)
{
    sql <- .myBindStatement(stmt, data)
    if (length(sql) > 1) {
        # parse out lists of values, collapse into one statement
        m <- sub('\\s+values\\s*(\\(.*\\))', '@\\1',
                 sql[-1], ignore.case=TRUE)
        m <- paste(sub('^[^@]+@', '', m), collapse=',')
        sql <- paste(sql[1], m, sep=',')
    }
    dbSendQuery(stmt$db, sql, ...)
}

.my.sql.exec <- function(db, sql, ..., chunk.kb=256, progress=FALSE)
{
    efn <- function(e)
    {
        .dbClearAll(db)
        stop(e)
    }
    data <- .sql.prep.data(...)
    sql <- .fixup.sql(db, sql)
    nr <- 0

    ps <- tryCatch(.myPrepareStatement(db, sql, data), error=efn)
    nd <- nrow(data)
    if (nd && progress) progress.bar(0, nd)
    if (nd && (regexpr('^insert.*\\svalues\\s*\\(', sql) > 0)) {
        chunk.rows <- .chunk.rows(data, chunk.kb)
        for (lo in seq(1,nd,chunk.rows)) {
            d <- data[lo:min(lo+chunk.rows-1,nd),,drop=FALSE]
            rs <- tryCatch(.myInsertMultiple(ps, d), error=efn)
            nr <- nr + dbGetRowsAffected(rs)
            if (progress) progress.bar(nr, nd)
        }
    } else {
        for (i in 1:max(1,nd)) {
            rs <- tryCatch(.myExecStatement(ps, data[i,,drop=FALSE]),
                           error=efn)
            nr <- nr + dbGetRowsAffected(rs)
            if (nd && progress) progress.bar(i, nd)
        }
    }
    nr
}

.lite.sql.exec <- function(db, sql, ..., chunk.kb=256, progress=FALSE)
{
    efn <- function(e)
    {
        dbRollback(db)
        .dbClearAll(db)
        stop(e)
    }
    sql <- .fixup.sql(db, sql)
    sql <- gsub(":[0-9]+", "?", sql)

    if (!length(list(...))) {
        dbBeginTransaction(db)
        rs <- tryCatch(dbSendQuery(db, sql), error=efn)
        nr <- dbGetRowsAffected(rs)
        dbCommit(db)
        return(nr)
    }

    data <- .sql.prep.data(...)
    nu <- nr <- 0
    nd <- nrow(data)
    if (!nd) return(0)
    chunk.rows <- .chunk.rows(data, chunk.kb)
    dbBeginTransaction(db)
    if (progress) progress.bar(0, nd)
    for (lo in seq(1,nd,chunk.rows)) {
        d <- data[lo:min(lo+chunk.rows-1,nd),,drop=FALSE]
        rs <- tryCatch(dbSendPreparedQuery(db, sql, d), error=efn)
        nr <- nr + dbGetRowsAffected(rs)
        nu <- nu + nrow(d)
        if (progress) progress.bar(nu, nd)
    }
    dbCommit(db)
    nr
}

sql.exec <- function(db, sql, ..., chunk.kb=256, progress=FALSE)
standardGeneric('sql.exec')

setGeneric('sql.exec', sql.exec)
setMethod('sql.exec', 'OraConnection', def=.ora.sql.exec)
setMethod('sql.exec', 'MySQLConnection', def=.my.sql.exec)
setMethod('sql.exec', 'SQLiteConnection', def=.lite.sql.exec)

