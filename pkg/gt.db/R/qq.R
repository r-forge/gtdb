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

prepanel.qqpval <-
function(x, n, groups=NULL, subscripts, ...)
{
    yy <- -log10(range(x, na.rm=TRUE)[c(2,1)])
    if (is.na(n)) {
        if (is.null(groups)) {
            n <- sum(!is.na(x))
        } else {
            sx <- split(x, groups[subscripts])
            n <- max(sapply(sx, function(x) sum(!is.na(x))))
        }
    }
    xx <- -log10(c(sum(!is.na(x))/n, 1/n))
    list(xlim=xx, ylim=yy, dx=diff(xx), dy=diff(yy))
}

panel.qqpval <-
function(x, n, max.pts, groups=NULL, ...)
{
    ffn <- function(n,z,a)
    {
        p <- ppoints(n,a)
        if (z==n) return(p)
        p <- rev(p[1:z])
        q <- 10^(seq(0,log10(p[1]),length.out=max.pts-z+1))
        1-c(q,p[-1])
    }
    dfn <- function(p) -log10((1-p)*length(x)/n)
    if (!is.null(groups)) {
        panel.superpose(x, y=NULL, n=n, max.pts=max.pts,
            groups=groups, panel.groups = panel.qqpval, ...)
    } else {
        if (is.na(n)) n <- sum(!is.na(x))
        x <- -log10(x[!is.na(x)])
        z <- which.min(sort(x-dfn(0),decr=TRUE)[1:max.pts] / (max.pts:1))
        yy <- quantile(x, ffn(length(x),z,1))
        xx <- dfn(ffn(length(x),z,0.5))
        panel.qq(xx, yy, ...)
    }
}

qqpval <-
function(x, ..., n=NA, max.pts=1000, prepanel=prepanel.qqpval,
         panel=panel.qqpval, xlab='theoretical quantiles')
{
    # this hackery is to give us a Q-Q plot with inverted log axes
    xs <- function(lim, logsc, ...)
    {
        s <- xscale.components.default(lim=-lim, logsc=TRUE, ...)
        s$num.limit <- -s$num.limit
        s$bottom$ticks$at <- -s$bottom$ticks$at
        s
    }
    ys <- function(lim, logsc, ...)
    {
        s <- yscale.components.default(lim=-lim, logsc=TRUE, ...)
        s$num.limit <- -s$num.limit
        s$left$ticks$at <- -s$left$ticks$at
        s
    }
    qqmath(x, ..., n=n, max.pts=max.pts, prepanel=prepanel, panel=panel,
           xscale.components=xs, yscale.components=ys, xlab=xlab)
}

panel.qqthin <-
function(x, max.pts, groups=NULL, ...)
{
    ffn <- function(n,z,a)
    {
        p <- ppoints(n,a)
        if (2*z >= n) return(p)
        p1 <- p[1:z]
        p2 <- p[(n-z+1):n]
        p <- ppoints(n,0.5)
        q <- pnorm(seq(qnorm(p[z]), qnorm(p[n-z+1]),
                       length.out=max.pts-2*z+2))
        c(p1,q[c(-1,-length(q))],p2)
    }
    if (!is.null(groups)) {
        panel.superpose(x, y=NULL, max.pts=max.pts,
                        groups=groups, panel.groups=panel.qqthin, ...)
    } else {
        x <- x[!is.na(x)]
        z <- which.min(sort(abs(x), decr=TRUE)[1:max.pts] / (max.pts:1))
        z <- ceiling(0.5*z)
        yy <- quantile(x, ffn(length(x),z,1))
        xx <- qnorm(ffn(length(x),z,0.5))
        ref <- trellis.par.get("reference.line")
        panel.qqmathline(x, col=ref$col, lty=ref$lty, lwd=ref$lwd)
        panel.xyplot(xx, yy, ...)
    }
}

qqthin <- function(x, ..., max.pts=1000, panel=panel.qqthin)
{
    qqmath(x, ..., panel=panel, max.pts=max.pts)
}
