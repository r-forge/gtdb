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

.flag.scope <-
  c('assay'='platform',
    'assay_position'='mapping',
    'assay_data'='dataset')

ls.flag <- function(target, parent.name)
{
    parent.id <- lookup.id(.flag.scope[target], parent.name)
    sql <-
     'select position, name flag_name, description, created_by, created_dt
      from %s_flag where %s_id=:1'
    sql <- sprintf(sql, target, .flag.scope[target])
    sql.query(gt.db::.gt.db, sql, parent.id)
}

mk.flag <-
function(target, parent.name, position, name, description)
{
    parent.id <- lookup.id(.flag.scope[target], parent.name)
    sql <-
     'insert into %s_flag
      values (:1,:2,:3,:4,:user:,:sysdate:)'
    sql <- sprintf(sql, target)
    sql.exec(gt.db::.gt.db, sql, parent.id, position, name, description)
}

#---------------------------------------------------------------------

unpack.flags <- function(target, parent.name, flags, ...)
{
    f <- ls.flag(target, parent.name, ...)
    if (nrow(f) == 0)
        return(data.frame(row.names=1:length(flags)))
    r <- sapply(2^(f$position-1),
                function(x) as.logical((flags %/% x) %% 2))
    structure(as.data.frame(r), names=f$flag.name)
}

pack.flags <- function(target, parent.name, data)
{
    f <- ls.flag(target, parent.name)
    data <- scale(data.matrix(data[,f$name]), 
                  center=FALSE, scale=2^(f$position-1))
    rowSums(data)
}
