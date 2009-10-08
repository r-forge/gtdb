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

# Miscellaneous support functions that don't fit elsewhere

#---------------------------------------------------------------------

.options <-
    list(db.mode='raw',
         tx.mode='hex')

.gt.db.options <- function(name, value, ...)
{
    opts <- gt.db::.options
    if (!missing(value)) {
        x <- list(value)
        names(x) <- name
    } else if (!missing(name) && is.list(name)) {
        x <- name
    } else if (!missing(name)) {
        return(opts[[name]])
    } else {
        x <- list(...)
        if (!length(x))
            return(opts)
    }
    opts[names(x)] <- x
    assign('.options', opts, 'package:gt.db')
    invisible()
}
