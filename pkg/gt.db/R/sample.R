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

# Functions for managing data about dataset samples

#---------------------------------------------------------------------

.fixup.gender <- function(gender)
{
    if (!all(gender %in% c('F','M',NA)))
        stop('Invalid gender data')
    factor(gender, levels=c('F','M'))
}

ls.sample <- function(dataset.name, show.ids=FALSE)
{
    sql <-
     'select sample_id, s.name sample_name, u.subject_id,
             u.name subject_name, gender, position
      from sample s, subject u
      where s.dataset_id=:1 and s.subject_id=u.subject_id'
    r <- sql.query(gt.db::.gt.db, sql, lookup.id('dataset', dataset.name))
    r$gender <- .fixup.gender(r$gender)
    .filter.ids(data.frame(r, row.names=r$sample.name), show.ids)
}

mk.sample <- function(dataset.name, data)
{
    .check.name(data$sample.name)
    sql <- 'insert into sample values (null,:1,:2,:3,:4,:5)'
    dset.id <- lookup.id('dataset', dataset.name)
    proj.id <- ls.dataset(, dataset.name, show.ids=TRUE)$project.id
    subj.id <- lookup.id('subject', data$subject.name, project.id=proj.id)
    if (is.null(data$gender)) {
        warning("gender information is missing")
        data$gender <- NA
    }
    data$gender <- .fixup.gender(data$gender)
    if (is.null(data$position)) {
        warning("sample positions are missing")
        data$position <- NA
    }
    if (any(data$position < 1) ||
        any(data$position != round(data$position))) {
        stop("sample positions should be 1-based integers")
    }
    sql.exec(gt.db::.gt.db, sql, dset.id, subj.id,
             data[c('sample.name','gender','position')])
}

rm.sample <- function(dataset.name, sample.name)
{
    dset.id <- lookup.id('dataset', dataset.name)
    samp.id <- lookup.id('sample', sample.name, dataset.id=dset.id)
    sql <- 'delete from sample where sample_id=:1'
    sql.exec(gt.db::.gt.db, sql, samp.id)
}

up.sample <- function(dataset.name, data)
{
    sql <-
     'update sample set subject_id=:1, gender=:2, position=:3
      where sample_id=:4'
    dset.id <- lookup.id('dataset', dataset.name)
    samp.id <- lookup.id('sample', data$sample.name, dataset.id=dset.id)
    proj.id <- ls.dataset(, dataset.name, show.ids=TRUE)$project.id
    subj.id <- lookup.id('subject', data$subject.name, project.id=proj.id)
    sql.exec(gt.db::.gt.db, sql, subj.id,
             data[c('gender','position')], samp.id)
}

mk.sample.attr <-
function(dataset.name, data, description=names(data))
{
    d <- .encode.attr(data)
    mk.attr('sample', dataset.name, names(data),
            d$datatype, d$levels, description)
}

rm.sample.attr <- function(dataset.name, attr.name)
rm.attr('sample', dataset.name, attr.name)

ls.sample.attr <- function(dataset.name)
ls.attr('sample', dataset.name)

#---------------------------------------------------------------------

fetch.sample.data <-
function(dataset.name, cols, show.all=FALSE)
{
    dataset.id <- lookup.id('dataset', dataset.name)
    sql <-
     "select s.name id, a.name time, value
      from sample s, sample_attr a, sample_value v
      where s.dataset_id = :1
        and s.sample_id=v.sample_id
        and a.sample_attr_id=v.sample_attr_id"
    if (!missing(cols)) {
        sql <- paste(sql, 'and a.name=:2')
        d <- lapply(cols, function (x)
                    sql.query(gt.db::.gt.db, sql, dataset.id, x))
        d <- do.call(rbind, d)
    } else {
        sql <- paste(sql, 'and a.is_hidden <= :2')
        d <- sql.query(gt.db::.gt.db, sql, dataset.id, show.all)
    }

    s <- ls.sample(dataset.name)
    if (nrow(d)) {
        d <- reshape(d, direction='wide')
        names(d) <- sub('^value.','',names(d))
        s <- cbind(s, d[match(s$sample.name,d$id),-1,drop=FALSE])
    }
    .decode.attr(s, ls.attr('sample', dataset.name))
}

store.sample.data <- function(dataset.name, data)
{
    dset.id <- lookup.id('dataset', dataset.name)
    samp.id <- lookup.id('sample', data$sample.name, dataset.id=dset.id)
    for (n in c('sample.name','subject.name','gender','position')) {
        data[n] <- NULL
    }
    attr.id <- lookup.id('sample_attr', names(data), dataset.id=dset.id)
    sql <- 'insert into sample_value values (:1,:2,:3)'
    insert.column <- function(attr.id, col)
    {
        w <- !is.na(col)
        sql.exec(gt.db::.gt.db, sql, samp.id[w], attr.id, col[w])
    }
    mapply(insert.column, attr.id, data)
}
