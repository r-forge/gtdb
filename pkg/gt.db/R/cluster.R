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

panel.cluster <-
function(x, y, group.number=1, bounds=c(), min.points=4, ...)
{
    require(cluster)
    panel.xyplot(x, y, ...)
    if (!length(bounds))
        return()
    grp <- function(x,n) lapply(x, function(y) y[(n-1)%%length(y)+1])
    sl <- grp(trellis.par.get('superpose.line'), group.number)
    if (!sl$lty)
        return()
    w <- !is.na(x) & !is.na(y)
    if (sum(w) < min.points)
        return()
    Cxy <- cov.wt(cbind(x,y)[w,])
    for (d2 in qchisq(bounds, df=2)) {
        llines(ellipsoidPoints(Cxy$cov, d2, Cxy$center),
               col=sl$col, alpha=sl$alpha, lty=sl$lty, lwd=sl$lwd)
    }
}

.gt.settings <- list(
    superpose.symbol=list(
        pch=c(1,1,1,3),
        col=c('#377db8','#e31a1c','#2daf4a','black')
    ),
    superpose.line=list(
        lty=c(1,1,1,0),
        col=c('#377db8','#e31a1c','#2daf4a','black')
    )
)

gt.cluster.plot <-
function(data, rescale=TRUE, bounds=c(0.5,0.95), min.points=4,
         between=list(x=0.5,y=0.5), scales=list(alternating=0),
         xlab=NULL, ylab=NULL, par.settings=.gt.settings, ...)
{
    if ('signal.a' %in% names(data)) {
        x <- data$signal.a
        y <- data$signal.b
        equal <- TRUE
    } else if ('fwd.a' %in% names(data)) {
        x <- data$fwd.a+data$rev.a
        y <- data$fwd.b+data$rev.b
        equal <- TRUE
    } else if ('strength' %in% names(data)) {
        x <- data$log.ratio
        y <- data$strength
        equal <- FALSE
    } else {
        stop("raw data not available")
    }

    if (equal) {
        # Equal scaling of X and Y
        prepanel <- function(...)
        {
            p <- prepanel.default.xyplot(...)
            p$xlim <- p$ylim <- c(min(p$xlim,p$ylim), max(p$xlim,p$ylim))
            p
        }
    } else {
        # center X at 0
        prepanel <- function(...)
        {
            p <- prepanel.default.xyplot(...)
            p$xlim <- c(-1,1)*max(abs(p$xlim))
            p
        }
    }

    n <- factor(data$assay.name)
    if (rescale) {
        if (equal) {
            q <- tapply(c(x,y), rep(n,2), range, na.rm=TRUE)
            q <- do.call('rbind', q)
            x <- 0.01 + 0.98*(x - q[n,1])/(q[n,2]-q[n,1])
            y <- 0.01 + 0.98*(y - q[n,1])/(q[n,2]-q[n,1])
        } else {
            q <- tapply(x, n, function(x) max(abs(x),na.rm=TRUE))
            x <- x / q[n]
            q <- tapply(y, n, range, na.rm=TRUE)
            q <- do.call('rbind', q)
            y <- 0.01 + 0.98*(y - q[n,1])/(q[n,2]-q[n,1])
        }
    }

    gt <- data$genotype
    if (is.numeric(gt)) {
        gt <- factor(gt, levels=0:3)
        gt[is.na(gt)] <- '3'
    }
    if (length(levels(gt)) == 3) {
        l <- levels(gt)
        nn <- gsub('[A-Z]','N',gsub('[a-z]','n', l[1]))
        gt <- factor(gt, levels=c(l, nn))
        gt[is.na(gt)] <- nn
    }
    p <- xyplot(y~x|n, groups=gt, bounds=bounds, min.points=min.points,
                scales=scales, prepanel=prepanel, aspect=1,
                panel=panel.superpose, panel.groups=panel.cluster,
                xlab=xlab, ylab=ylab, between=between,
                par.settings=par.settings, ...)
    p
}

xyplot.gt.data <- function(x, ...)
{
    gt.cluster.plot(reshape.gt.data(x, na.codes='n'), ...)
}

.multi.identify <-
function (radius = 6, panel.args = trellis.panelArgs())
{
    xy <- xy.coords(panel.args$x, panel.args$y)
    px <- convertX(unit(xy$x, "native"), "points", TRUE)
    py <- convertY(unit(xy$y, "native"), "points", TRUE)
    unmarked <- rep(TRUE, length(xy$x))
    while (sum(unmarked)) {
        ll <- grid.locator(unit = "points")
        if (is.null(ll))
            break
        lx <- convertX(ll$x, "points", TRUE)
        ly <- convertY(ll$y, "points", TRUE)
        pdists <- sqrt((px - lx)^2 + (py - ly)^2)
        w <- which(unmarked & (pdists < radius))
        grid.circle(lx, ly, radius, 'points',
                    gp=gpar(col='#a0a0a0', fill='#a0a0a0'))
        unmarked[w] <- FALSE
    }
    which(!unmarked)
}

adjust.gt.calls <-
function(data, ..., radius = 6)
{
    if (!is.factor(data$genotype))
        stop("genotype column should be a factor")
    if (length(unique(data$assay.name))!=1)
        stop('one SNP at a time please')
    if (is.null(data$orig.genotype))
        data$orig.genotype <- data$genotype
    opt <- ''
    levs <- paste(levels(data$genotype), collapse="','")
    if (any(is.na(data$genotype))) levs <- c(levs, 'NA')
    prompt <- paste("Enter genotype to assign ('", levs,
                    "') or 'q' to exit: ", sep='')

    while (opt != 'q') {
        if (opt %in% levels(data$genotype)) {
            trellis.focus("panel", 1, 1)
            cat("Select points to assign genotypes, then right click 'Stop'\n")
            data$genotype[.multi.identify(radius)] <- na.if(opt,'NA')
        } else {
            cat("Please enter a valid level.\n")
        }
        print(gt.cluster.plot(data, ...))
        print(table(data$genotype))
        cat("\n")
        opt <- readline(prompt)
    }
    data
}
