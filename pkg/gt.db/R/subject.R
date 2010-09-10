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

# Functions for managing data about project subjects

#---------------------------------------------------------------------

ls.subject <- function(project.name, show.ids=FALSE)
{
    sql <-
     'select subject_id, name subject_name
      from subject where project_id=:1'
    r <- sql.query(gt.db::.gt.db, sql, lookup.id('project', project.name))
    .filter.ids(data.frame(r, row.names=r$subject.name), show.ids)
}

mk.subject <- function(project.name, data)
{
    .check.name(data$subject.name)
    sql <- 'insert into subject values (null,:1,:2)'
    proj.id <- lookup.id('project', project.name)
    sql.exec(gt.db::.gt.db, sql, proj.id, data$subject.name)
}

rm.subject <- function(project.name, subject.name)
{
    proj.id <- lookup.id('project', project.name)
    subj.id <- lookup.id('subject', subject.name, project.id=proj.id)
    sql <- 'delete from subject where subject_id=:1'
    sql.exec(gt.db::.gt.db, sql, subj.id)
}

mk.subject.attr <-
function(project.name, data, description=names(data))
{
    d <- .encode.attr(data)
    mk.attr('subject', project.name, names(data),
            d$datatype, d$levels, description)
}

rm.subject.attr <- function(project.name, attr.name)
rm.attr('subject', project.name, attr.name)

ls.subject.attr <- function(project.name)
ls.attr('subject', project.name)

#---------------------------------------------------------------------

fetch.subject.data <-
function(project.name, cols, show.all=FALSE)
{
    project.id <- lookup.id('project', project.name)
    sql <-
     "select s.name id, a.name time, value
      from subject s, subject_attr a, subject_value v
      where s.project_id = :1
        and s.subject_id=v.subject_id
        and v.subject_attr_id=a.subject_attr_id"
    if (!missing(cols)) {
        sql <- paste(sql, 'and a.name=:2')
        d <- lapply(cols, function (x)
                    sql.query(gt.db::.gt.db, sql, project.id, x))
        d <- do.call(rbind, d)
    } else {
        sql <- paste(sql, 'and a.is_hidden <= :2')
        d <- sql.query(gt.db::.gt.db, sql, project.id, show.all)
    }

    s <- ls.subject(project.name)
    if (nrow(d)) {
        d <- reshape(d, direction='wide')
        names(d) <- sub('^value.','',names(d))
        s <- cbind(s, d[match(s$subject.name,d$id),-1,drop=FALSE])
    }
    .decode.attr(s, ls.attr('subject', project.name))
}

store.subject.data <- function(project.name, data)
{
    proj.id <- lookup.id('project', project.name)
    sn <- data$subject.name ; data$subject.name <- NULL
    subj.id <- lookup.id('subject', sn, project.id=proj.id)
    attr.id <- lookup.id('subject_attr', names(data), project.id=proj.id)
    sql <- 'insert into subject_value values (:1,:2,:3)'
    insert.column <- function(attr.id, col)
    {
        w <- !is.na(col)
        if (!any(w)) return(0)
        sql.exec(gt.db::.gt.db, sql, subj.id[w], attr.id, col[w])
    }
    mapply(insert.column, attr.id, data)
}
