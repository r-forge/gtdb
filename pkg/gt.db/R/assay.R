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

ls.mapping <-
function(platform.name, show.all=FALSE, show.ids=FALSE)
{
    sql <-
     'select mapping_id, name mapping_name, description,
        is_hidden, created_by, created_dt
      from mapping
      where platform_id=:1
        and is_hidden<=:2'
    r <- sql.query(gt.db::.gt.db, sql, lookup.id('platform', platform.name),
                   show.all)
    .filter.ids(r, show.ids)
}

mk.mapping <-
function(platform.name, mapping.name, description, is.hidden=FALSE)
{
    .check.name(mapping.name)
    plat.id <- lookup.id('platform', platform.name)
    sql <-
     "insert into mapping
      values (null, :1, :2, :3, :4, :user:, :sysdate:)"
    sql.exec(gt.db::.gt.db, sql, plat.id, mapping.name,
             description, is.hidden)
}

rm.mapping <- function(platform.name, mapping.name)
{
    grp.id <- lookup.id('mapping', mapping.name,
                        platform.id=lookup.id('platform', platform.name))
    sql <- 'delete from mapping where mapping_id=:1'
    sql.exec(gt.db::.gt.db, sql, grp.id)
}

lookup.mapping.id <- function(platform.name, mapping.name)
{
    if (missing(mapping.name)) {
        m <- ls.mapping(platform.name, show.ids=TRUE)
        if (nrow(m) != 1)
            stop('could not choose default mapping', call.=FALSE)
        structure(m$mapping.id, names=m$mapping.name)
    } else {
        plat.id <- lookup.id('platform', platform.name)
        lookup.id('mapping', mapping.name, platform.id=plat.id)
    }
}

#---------------------------------------------------------------------

ls.assay.group <- function(platform.name, show.ids=FALSE)
{
    sql <-
     'select assay_group_id, name assay_group,
        description, created_by, created_dt
      from assay_group where platform_id=:1'
    r <- sql.query(gt.db::.gt.db, sql, lookup.id('platform', platform.name))
    .filter.ids(r, show.ids)
}

mk.assay.group <-
function(platform.name, assay.group, description)
{
    .check.name(assay.group)
    plat.id <- lookup.id('platform', platform.name)
    sql <-
     "insert into assay_group
      values (null, :1, :2, :3, :user:, :sysdate:)"
    sql.exec(gt.db::.gt.db, sql, plat.id, assay.group, description)
}

rm.assay.group <- function(platform.name, assay.group)
{
    grp.id <- lookup.id('assay_group', assay.group,
                        platform.id=lookup.id('platform', platform.name))
    sql <- 'delete from assay_group where assay_group_id=:1'
    sql.exec(gt.db::.gt.db, sql, grp.id)
}

#---------------------------------------------------------------------

ls.assay <- function(platform.name, assay.group='%', show.ids=FALSE)
{
    sql <-
     'select g.assay_group_id, g.name assay_group,
        assay_id, a.name assay_name, alleles, probe_seq
      from assay_group g, assay a
      where g.assay_group_id=a.assay_group_id
        and g.platform_id=:1
        and g.name like :2'
    r <- sql.query(gt.db::.gt.db, sql, lookup.id('platform', platform.name),
                   assay.group)
    .filter.ids(data.frame(r, row.names=r$assay.name), show.ids)
}

mk.assay <- function(platform.name, data)
{
    plat.id <- lookup.id('platform', platform.name)
    grp.id <- lookup.id('assay_group', unique(data$assay.group),
                        platform.id=plat.id)
    .check.name(data$assay.name)
    sql <- 'insert into assay values (null,:1,:2,:3,:4,:5)'
    sql.exec(gt.db::.gt.db, sql, plat.id, grp.id[data$assay.group],
             data[c('assay.name','alleles','probe.seq')])
}

#---------------------------------------------------------------------

.fixup.ploidy <- function(ploidy)
{
    if (!all(ploidy %in% c('A','M','X','Y',NA)))
        stop('Invalid ploidy data')
    factor(ploidy, levels=c('A','M','X','Y'))
}

ls.assay.position <-
function(platform.name, mapping.name, show.ids=FALSE)
{
    map.id <- lookup.mapping.id(platform.name, mapping.name)
    sql <-
     'select a.assay_id, a.name assay_name, scaffold, position,
        strand, ploidy, dbsnp_rsid, dbsnp_orient
      from assay a, assay_position p
      where mapping_id=:1
        and a.assay_id=p.assay_id'
    r <- sql.query(gt.db::.gt.db, sql, map.id)
    r$ploidy <- .fixup.ploidy(r$ploidy)
    .filter.ids(data.frame(r,row.names=r$assay.name), show.ids)
}

mk.assay.position <- function(platform.name, mapping.name, data)
{
    map.id <- lookup.mapping.id(platform.name, mapping.name)
    if (is.null(data$assay.id)) {
        plat.id <- lookup.id('platform', platform.name)
        data$assay.id <- lookup.id('assay', data$assay.name,
                                   platform.id=plat.id)
    }
    sql <-
     'insert into assay_position values (:1,:2,:3,:4,:5,:6,:7,:8)'
    if (is.null(data$ploidy)) {
        warning("ploidy information is missing")
        data$ploidy <- NA
    }
    data$ploidy <- .fixup.ploidy(data$ploidy)
    sql.exec(gt.db::.gt.db, sql, map.id,
             data[c('assay.id', 'scaffold', 'position', 'strand',
                    'ploidy', 'dbsnp.rsid', 'dbsnp.orient')])
}

#---------------------------------------------------------------------

mk.assay.data <- function(dataset.name, data)
{
    dset.id <- lookup.id('dataset', dataset.name)
    if (is.null(data$flags)) data$flags <- 0
    if (is.null(data$qscore)) data$qscore <- NA
    if (is.null(data$raw.data)) data$raw.data <- NA
    if (is.null(data$assay.id)) {
        plat.id <- ls.dataset(dataset.name=dataset.name,
                              show.ids=TRUE)$platform.id
        data$assay.id <- lookup.id('assay', data$assay.name,
                                   platform.id=plat.id)
    }

    db.mode <- .gt.db.options('db.mode')
    tx.mode <- .gt.db.options('tx.mode')
    if (db.mode == tx.mode) {
        cvt.fn <- ''
    } else if (db.mode == 'raw' && tx.mode == 'hex') {
        cvt.fn <- ':unhex:'
    } else {
        stop('unknown conversion!')
    }

    sql <-
     'insert into assay_data
      values (null,:1,:2,:3,:4,%1$s(:5),%1$s(:6))'
    sql <- sprintf(sql, cvt.fn)
    cols <- c('assay.id','flags','genotype','qscore','raw.data')
    sql.exec(gt.db::.gt.db, sql, dset.id, data[cols])
}
