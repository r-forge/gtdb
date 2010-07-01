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

gplot <-
function(formula, data, aggr.fn=max, rescale=FALSE, binsz=1e6,
         subset=TRUE, col.regions=rev(heat.colors(100)[10:90]),
         scales=list(x=list(at=seq(0,250,20), draw=TRUE)),
         shrink=list(x=1,y=0.75), colorkey=list(height=0.25),
         xlab='Position, Mb', ylab='Chromosome', zlim, ...)
{
    subset <- eval(substitute(subset), data, environment(formula))
    data <- data[subset,]
    z <- latticeParseFormula(formula, data)$right
    if (!missing(zlim))
        z <- pmax(zlim[1],pmin(zlim[2],z))
    bin <- list(factor(data$scaffold), trunc(data$position/binsz))
    m <- tapply(z, bin, aggr.fn)
    if (rescale) m <- m/mean(m,na.rm=TRUE)
    colnames(m) <- (as.numeric(colnames(m))+0.5) * binsz / 1e6

    chr <- .sort.levels(rownames(m), decreasing=TRUE)
    pos <- rep(as.numeric(colnames(m)), each=nrow(m))
    d <- data.frame(pos=pos, chr=chr, z=as.vector(m))
    p <- function(...) {
        panel.abline(v=scales$x$at,col='gray');
        panel.levelplot(...)
    }

    levelplot(z~pos+chr, d, scales=scales, aspect='fill',
              col.regions=col.regions, colorkey=colorkey,
              shrink=shrink, xlab=xlab, ylab=ylab, panel=p, ...)
}

manhattan.plot <-
    function(y, data, gap=0, threshold=-log10(5e-8), around=0,
             xticks=c(1:12,14,16,18,20,22,'X','Y'), cex=0.25,
             xlab=NULL, ylab=deparse(substitute(y)),
             col=c('#d0d0d0','#e0e0e0','#ff0000'), ...)
{
    val <- eval(substitute(y), data, parent.frame())
    len <- with(data, tapply(position, .sort.levels(scaffold),
                             max, na.rm=TRUE))
    ofs <- cumsum(c(0,as.numeric(len+gap)))
    mid <- (ofs[-1] + ofs[-length(ofs)])/2
    pos <- data$position + ofs[match(data$scaffold, names(len))]
    set <- list(superpose.symbol=list(col=col))
    grp <- (match(data$scaffold, names(len)) %% 2)
    grp <- ifelse(val > threshold, 2, grp)
    if (around > 0) {
        hits <- pos[which(val > threshold)]
        keep <- (diff(hits) > around)
        hits <- hits[c(TRUE,keep) | c(keep,TRUE)]
        for (hit in hits) {
            grp[abs(pos - hit) < around] <- 2
        }
    }
    xlim <- range(pos,na.rm=TRUE)
    xlim <- xlim + (xlim[2]-xlim[1]) * c(-0.02,0.02)
    panel.fn <- function(...) { xyplot(...) ; panel.refline(h=threshold) }
    if (length(len) == 1) {
        set <- list(superpose.symbol=list(col=col[-2]))
        xyplot(val~pos, ..., cex=cex, groups=grp, par.settings=set,
               panel=panel.fn, xlab=xlab, ylab=ylab, xlim=xlim)
    } else {
        w <- match(paste('chr',xticks,sep=''), names(len))
        xyplot(val~pos, ..., cex=cex, groups=grp, par.settings=set,
               scales=list(x=list(at=mid[w],tck=0,labels=xticks)),
               panel=panel.fn, xlab=xlab, ylab=ylab, xlim=xlim)
    }
}
