#
# Copyright (C) 2009, Perlegen Sciences, Inc.
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

# functions for managing the attribute tables

.attr.scope <- c(subject='project', sample='dataset')

ls.attr <-
function(target, parent.name, show.all=FALSE, show.ids=FALSE)
{
    sql <-
     'select %1$s_attr_id, name attr_name, datatype, levels,
        description, is_hidden, created_by, created_dt
      from %1$s_attr
      where %2$s_id = :1 and is_hidden <= :2'
    sql <- sprintf(sql, target, .attr.scope[target])
    parent.id <- lookup.id(.attr.scope[target], parent.name)
    r <- sql.query(gt.db::.gt.db, sql, parent.id, show.all)
    .filter.ids(data.frame(r, row.names=r$attr.name), show.ids)
}

mk.attr <-
function(target, parent.name, attr.name, datatype,
         levels, description, is.hidden=FALSE)
{
    parent.id <- lookup.id(.attr.scope[target], parent.name)
    .check.name(attr.name)
    sql <-
     'insert into %s_attr
      values (null,:1,:2,:3,:4,:5,:6,:user:,:sysdate:)'
    sql <- sprintf(sql, target)
    sql.exec(gt.db::.gt.db, sql, parent.id, attr.name, datatype,
             levels, description, is.hidden)
}

rm.attr <- function(target, parent.name, attr.name)
{
    sql <- 'delete from %s_attr where %s_id=:1 and name=:2'
    sql <- sprintf(sql, target, .attr.scope[target])
    parent.id <- lookup.id(.attr.scope[target], parent.name)
    sql.exec(gt.db::.gt.db, sql, parent.id, attr.name)
}

#---------------------------------------------------------------------

# Functions for converting attribute data between R and database

.decode.attr <- function(data, info)
{
    for (n in intersect(names(data), info$attr.name)) {
        dt <- info[n,'datatype']
        if (dt == 'number') {
            data[,n] <- as.numeric(data[,n])
        } else if (dt == 'boolean') {
            data[,n] <- as.logical(as.numeric(data[,n]))
        } else if (dt == 'factor') {
            lv <- sub("^'(.*)'$", "\\1", info[n,"levels"])
            lv <- strsplit(lv, "','")[[1]]
            data[,n] <- factor(data[,n], levels=lv)
        }
    }
    data
}

.encode.attr <- function(data, drop=c())
{
    data <- data[!(names(data) %in% drop)]
    fn <- function(x)
        if (is.factor(x)) {
            lv <- paste(levels(x),collapse="','")
            c('factor', paste("'", lv, "'", sep=''))
        } else if (is.logical(x)) {
            c('boolean', "'0','1'")
        } else {
            c(if (is.numeric(x)) 'number' else 'string', NA)
        }
    r <- do.call('rbind', lapply(data, fn))
    data.frame(attr.name=names(data), datatype=r[,1], levels=r[,2],
               row.names=names(data), stringsAsFactors=FALSE)
}

