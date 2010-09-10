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

# these functions implement a subclass of data frames for which
# attributes are "sticky" and survive indexing

kept.attr <- function(x)
{
    a <- attributes(x)
    a[c('names','row.names','class','dim','dimnames')] <- NULL
    a
}

keep.attr <- function(.Data, ..., .Attr=NULL)
{
    cl <- union('keep.attr', class(.Data))
    do.call('structure', c(list(.Data, class=cl, ...), .Attr))
}

'[.keep.attr' <- function(x, ...)
{
    d <- NextMethod()
    if (identical(class(x), class(d)))
      keep.attr(d, .Attr=kept.attr(x))
    else
      d
}

'[<-.keep.attr' <- function(x, ..., value)
keep.attr(NextMethod(), .Attr=kept.attr(x))
