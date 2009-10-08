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

# functions for manipulating packed strings of T/F values

as.mask <- function(x)
{
    if (!is.character(x))
        x <- paste(ifelse(if.na(x,FALSE),'T','F'),collapse='')
    x
}

un.mask <- function(s)
{
    if (is.character(s))
        s <- as.logical(strsplit(s,'')[[1]])
    s
}

mask.str <- function(str, mask, ch='x')
{
    .Call('do_mask_str', str, as.mask(mask), ch, PACKAGE='gt.db')
}

# fast character-counting functions

nsubstr <- function(a, b) .Call('do_nsubstr', a, b, PACKAGE='gt.db')

ch.table <- function(s1, s2, chars)
{
    if (missing(s2))
        .Call('do_ch_table', s1, chars, PACKAGE='gt.db')
    else
        .Call('do_ch_table_2', s1, s2, chars, PACKAGE='gt.db')
}
