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

# Administrative functions for viewing and modifying the contents
# of the project, platform, and dataset tables

#---------------------------------------------------------------------

ls.project <-
function(project.name='%', show.all=FALSE, show.ids=FALSE)
{
    sql <-
     "select project_id, name project_name, description,
             (select count(*) from dataset d
              where d.project_id=p.project_id
                and is_hidden <= :1) datasets,
             is_hidden, created_by, created_dt
      from project p
      where name like :2
        and is_hidden <= :3"
    r <- sql.query(gt.db::.gt.db, sql, show.all, project.name, show.all)
    .filter.ids(r, show.ids)
}

mk.project <- function(project.name, description, is.hidden=FALSE)
{
    .check.name(project.name)
    sql <-
     "insert into project
      values (null, :1, :2, :3, :user:, :sysdate:)"
    sql.exec(gt.db::.gt.db, sql, project.name, description, is.hidden)
}

rm.project <- function(project.name)
{
    sql <- 'delete from project where project_id=:1'
    sql.exec(gt.db::.gt.db, sql, lookup.id('project', project.name))
}

#---------------------------------------------------------------------

ls.platform <- function(platform.name='%', show.ids=FALSE)
{
    sql <-
     "select platform_id, name platform_name, description,
             (select count(*) from assay_group g
              where p.platform_id=g.platform_id) assay_groups,
             (select count(*) from dataset d
              where p.platform_id=d.platform_id) datasets,
             created_by, created_dt
      from platform p
      where name like :1"
    r <- sql.query(gt.db::.gt.db, sql, platform.name)
    .filter.ids(r, show.ids)
}

mk.platform <- function(platform.name, description)
{
    .check.name(platform.name)
    sql <-
     "insert into platform
      values (null, :1, :2, :user:, :sysdate:)"
    sql.exec(gt.db::.gt.db, sql, platform.name, description)
}

rm.platform <- function(platform.name)
{
    sql <- 'delete from platform where platform_id=:1'
    sql.exec(gt.db::.gt.db, sql, lookup.id('platform', platform.name))
}

#---------------------------------------------------------------------

ls.dataset <-
function(project.name='%', dataset.name='%',
         show.all=FALSE, show.ids=FALSE)
{
    sql <-
     "select dataset_id, d.name dataset_name,
             d.project_id, p.name project_name,
             d.platform_id, m.name platform_name,
             d.description, d.raw_layout, d.is_hidden,
             d.created_by, d.created_dt
      from dataset d, project p, platform m
      where d.project_id=p.project_id
        and d.platform_id=m.platform_id
        and p.name like :1
        and d.name like :2
        and p.is_hidden <= :3
        and d.is_hidden <= :4"
    r <- sql.query(gt.db::.gt.db, sql, project.name, dataset.name,
                   show.all, show.all)
    .filter.ids(r, show.ids)
}

mk.dataset <-
function(dataset.name, project.name, platform.name,
         description, raw.layout=c(NA,'signal','seqread'),
         is.hidden=FALSE)
{
    raw.layout <- match.arg(raw.layout)
    .check.name(dataset.name)
    proj.id <- lookup.id('project', project.name)
    plat.id <- lookup.id('platform', platform.name)
    sql <-
     "insert into dataset
      values (null, :1, :2, :3, :4, :5, :6, :user:, :sysdate:)"
    sql.exec(gt.db::.gt.db, sql, proj.id, plat.id, dataset.name,
             description, raw.layout, is.hidden)
}

rm.dataset <- function(dataset.name)
{
    sql <- 'delete from dataset where dataset_id=:1'
    sql.exec(gt.db::.gt.db, sql, lookup.id('dataset', dataset.name))
}
