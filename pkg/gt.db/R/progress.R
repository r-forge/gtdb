#
# Copyright (C) 2009, Perlegen Sciences, Inc.
# Copyright (C) 2010, 23andMe, Inc.
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

progress.bar <- function (done, total, width=getOption('width'))
{
    len <- width - 40
    now <- Sys.time()
    p <- try(get(".Progress", pos='package:gt.db'), silent=TRUE)
    if ((class(p) == 'try-error') || (done == 0)) {
        bg <- paste(c("[", rep("-", len), "]"), collapse="")
        fg <- paste(c("[", rep("0", len), "]"), collapse="")
        p <- list(bg=bg,fg=fg,start=now)
    } else {
        if ((difftime(now,p$now) < 1) && (done < total)) {
            assign(".Progress", p, pos='package:gt.db')
            return(invisible())
        }
    }
    p$now <- now
    dt <- unclass(difftime(now,p$start,units='secs'))
    brk <- floor(len * done/total + 1.00001)
    t0 <- as.POSIXct('2001/01/01')
    t1 <- t0 + dt
    t2 <- t0 + total*dt/done - dt
    if (done == total) {
        cat(p$fg, " ", format(t1,format="%H:%M:%S"),
            " elapsed                     \n", sep="")
        rm(".Progress", pos='package:gt.db')
    } else {
        cat(substr(p$fg,1,brk), substr(p$bg,brk+1,len+2), " ",
            format(t1,format="%H:%M:%S"), " elapsed, ",
            format(t2,format="%H:%M:%S"), " remaining\r", sep="")
        assign(".Progress", p, pos='package:gt.db')
    }
    invisible(flush.console())
}
