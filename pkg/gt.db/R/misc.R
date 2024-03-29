#
# Copyright (C) 2009, Perlegen Sciences, Inc.
# Copyright (C) 2010, 23andMe, Inc.
#
# Written by David A. Hinds <dhinds@23andMe.com>
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

# Miscellaneous support functions that don't fit elsewhere

#---------------------------------------------------------------------

if.na <- function(val,yes,no=val) ifelse(is.na(val),yes,no)

na.if <- function(v1,v2) ifelse(v1==v2,NA,v1)

rawToHex <- function(raw)
.Call("raw_to_hex", as.vector(raw), PACKAGE="gt.db")

hexToRaw <- function(hex)
as.vector(.Call("hex_to_raw", hex[1], PACKAGE="gt.db"))

ramToHex <- function(raw)
.Call("raw_to_hex", raw, PACKAGE="gt.db")

hexToRam <- function(hex)
.Call("hex_to_raw", hex, PACKAGE="gt.db")

ramToChar <- function(raw)
.Call("raw_to_char", raw, PACKAGE="gt.db")

charToRam <- function(str)
.Call("char_to_raw", str, PACKAGE="gt.db")

#---------------------------------------------------------------------

# validation rules for names of objects

.check.name <- function(name)
{
    bad <- name[regexpr('^[A-Za-z0-9_.]+$', name) < 0]
    if (length(bad) == 1) {
        stop("invalid name '", bad, "'", call.=FALSE)
    } else if (length(bad)) {
        stop("invalid names '", bad[1], "', '", bad[2],
             ifelse(length(bad)>2, "', ...", "'"), call.=FALSE)
    }
}

#---------------------------------------------------------------------

revcomp <- function(x, ambig=FALSE)
{
    if (ambig)
        x <- chartr('ACGTMRWSYKVHDBacgtmrwsykvhdb',
                    'TGCAKYWSRMBDHVtgcakywsrmbdhv', x)
    else
        x <- chartr('ACGTacgt','TGCAtgca', x)
    sapply(x, function(s) rawToChar(rev(charToRaw(s))), USE.NAMES=FALSE)
}

#---------------------------------------------------------------------

use.gt.db <- function(dbConnection)
{
    assign('.gt.db', dbConnection, 'package:gt.db')
    x <- try(sql.query(gt.db::.gt.db, 'select * from gtdb_option'),
             silent=TRUE)
    if (class(x) == 'data.frame') {
        mapply(.gt.db.options, x$name, x$value)
    }
    invisible()
}

init.gt.db <- function(db.mode=c('raw','hex','zip'))
{
    db.mode <- match.arg(db.mode)
    path <- library(help='gt.db')$path
    schema <- switch(class(gt.db::.gt.db),
                     SQLiteConnection='mk_sqlite.sql',
                     MySQLConnection='mk_mysql.sql',
                     OraConnection='mk_oracle.sql')
    st <- readLines(paste(path, 'schema', schema, sep='/'))
    st <- st[!grepl('^--',st)]
    # split into statements at ';' followed by empty line
    st <- strsplit(paste(st, collapse='\n'), ';\n\n')[[1]]
    st <- st[!grepl('^\\s*delimiter\\s+',st)]
    sapply(st, sql.exec, db=gt.db::.gt.db, USE.NAMES=FALSE)
    .gt.db.options(db.mode=db.mode)
    sql.exec(gt.db::.gt.db, 'insert into gtdb_option values (:1,:2)',
             'db.mode', db.mode)
}

gt.demo.check <- function()
{
    prompt <- 'press <return> to create temporary demo database: '
    if (!exists('.gt.db','package:gt.db')) {
        e <- Sys.getenv('GT_DB_DEMO')
        if (e != '') {
            s <- strsplit(e, ',')[[1]]
            if (s[1] == 'MySQL') {
                require(RMySQL)
                db <- dbConnect(dbDriver('MySQL'), s[2], s[3], s[4])
            } else if (s[1] == 'SQLite') {
                require(RSQLite)
                db <- dbConnect(dbDriver('SQLite'), s[2],
                                loadable.extensions=TRUE)
                if (!is.na(s[3]))
                    sql <- sprintf("select load_extension('%s')", s[3])
                sql.query(db, sql)
            } else if (s[1] == 'Oracle') {
                require(ROracle)
                db <- dbConnect(dbDriver('Oracle'), s[2], s[3], s[4])
            }
            use.gt.db(db)
        } else if (!interactive() || !is.null(readline(prompt))) {
            require(RSQLite)
            warning("creating in-memory demo database...")
            db <- dbConnect(dbDriver('SQLite'), ':memory:',
                            loadable.extensions=TRUE)
            use.gt.db(db)
            init.gt.db(db.mode='hex')
            demo('setup.gt.demo')
        } else {
            stop('No GT.DB database connection', call.=FALSE)
        }
    }
    if (!nrow(ls.project('Demo')))
        stop('Demo project data is unavailable', call.=FALSE)
}

#---------------------------------------------------------------------

lookup.id <- function(table, name, ..., use.index)
{
    sql <- sprintf('select name, %1$s_id id from %1$s', table)
    args <- list(...)
    na <- length(args)
    where <- function(bind)
    { nm <- chartr('.','_',names(args))
      paste(nm, ' = :', bind, collapse=' and ', sep='') }
    if (missing(use.index))
        use.index <- (length(name) < 10)
    if (use.index) {
        sql <- paste(sql, 'where name=:1')
        if (na) sql <- paste(sql, 'and', where(1+(1:na)))
        find.one <- function(x, ...) sql.query(gt.db::.gt.db, sql, x, ...)
        x <- lapply(name, find.one, ..., USE.NAMES=FALSE)
        x <- do.call('rbind', x)
    } else {
        if (na) sql <- paste(sql, 'where', where(1:na))
        x <- sql.query(gt.db::.gt.db, sql, ...)
    }
    if (nrow(x)) {
        if (anyDuplicated(x[,1]))
            stop('lookup failed: please be more specific', call.=FALSE)
        r <- x[match(name,x[,1]),2]
        f <- name[is.na(r)]
    } else f <- name
    if (length(f) > 0)
        stop('lookup failed for ', table, " '", f[1], "'", call.=FALSE)
    structure(r,names=name)
}

.filter.ids <- function(data, show.ids, keep=c())
{
    if (show.ids) return(data)
    n <- names(data)
    data[(regexpr('\\.id$', n) < 0) | (n %in% keep)]
}

set.hidden <- function(table, name, is.hidden=TRUE, ...)
{
    id <- lookup.id(table, name, ...)
    sql <- 'update %i$s set is_hidden=:1 where %1$s_id=:2'
    sql.exec(gt.db::.gt.db, sprintf(sql, table), is.hidden, id)
}

#---------------------------------------------------------------------

fetch.pt.data <-
function(dataset.name, cols, pca=FALSE, show.all=FALSE)
{
    proj <- ls.dataset(, dataset.name)$project.name
    sa <- fetch.sample.data(dataset.name, cols, show.all)
    su <- fetch.subject.data(proj, cols, show.all)
    pt <- cbind(sa, su[sa$subject.name,-1,drop=FALSE])
    if (!identical(pca,FALSE)) {
        args <- c(list(dataset.name), if (is.list(pca)) pca)
        pt <- cbind(pt, do.call('fetch.prcomp', args)$loadings)
    }
    pt <- subset(pt, !is.na(position))
    keep.attr(pt[order(pt$position),], dataset.name=dataset.name)
}
