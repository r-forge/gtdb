#
# Copyright (C) 2010, 23andMe, Inc.
#
# Written by David A. Hinds <dhinds@23andMe.com>
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

setup.tracks <- function(xlim, name='region', newpage=TRUE, ...)
{
    if (newpage)
        grid.newpage()
    pushViewport(viewport(name=name, xscale=xlim, ...))
}

draw.tracks <- function(tracks, scale=c('Mb','Kb','bp'), xlab)
{
    tck <- function(x) { p <- pretty(x) ; p[p >= x[1] & p <= x[2]] }
    xscale <- current.viewport()$xscale
    scale <- match.arg(scale)
    mul <- switch(match.arg(scale), Mb=1e6, Kb=1e3, 1)
    grid.xaxis(tck(xscale), format(tck(xscale/mul)))
    if (missing(xlab))
        xlab <- paste(gt$scaffold[1], 'position,', scale)
    if (!is.null(xlab))
        grid.text(xlab, y=unit(-3,'lines'))
    y <- unit(0,'npc')
    for (tn in 1:length(tracks)) {
        name <- names(tracks)[tn]
        if (is.null(name) || name=='')
            name <- sprintf("track%d", tn)
        tr <- tracks[[tn]]
        args <- list(name=name, y=y, just='bottom', xscale=xscale)
        vp <- do.call('viewport', c(args, tr$vp))
        pushViewport(vp)
        if (tn == length(tracks))
            grid.xaxis(tck(xscale), FALSE, FALSE)
        grid.rect()
        grid.draw(tr$grob)
        upViewport()
        y <- y + tr$vp$height
    }
    pushViewport(viewport(name='border', y=unit(0,'npc'),
                          height=y, just='bottom'))
    grid.rect()
    upViewport()
}

ld.track <-
    function(gt, col=gray(seq(1,0,-0.01)), margin=0.02, fill='#f0f0f0')
{
    ld <- ld.gt.data(gt, measure='rsqr')
    dimnames(ld) <- list(gt$position,gt$position)
    xd <- as.data.frame.table(ld, stringsAsFactors=FALSE)
    xd[,1] <- as.numeric(xd[,1])
    xd[,2] <- as.numeric(xd[,2])
    names(xd) <- c('x','y','z')
    cv <- current.viewport()
    if (identical(cv$xscale,c(0,1)))
        stop('track coordinates not initialized.  Use setup.tracks()')
    yscale <- cv$xscale + margin*c(-1,1)*diff(cv$xscale)
    height <- cv$width * (0.5+margin)
    xd <- subset(xd, x>y)
    g <- .panel.ldplot(xd$x, xd$y, xd$z, at=seq(0,1,0.01),
                       col.regions=col, draw=FALSE)
    f <- rectGrob(gp=gpar(fill=fill))
    invisible(list(grob=grobTree(f,g),
                   vp=list(yscale=yscale, height=height)))
}

geneGrob <-
    function(name, strand, txStart, txEnd, col, ofs=0,
             arrow.thick=0.4, margin=0.2, text.sep=0.1)
{
    arrow <- function(x1, x2)
    {
        flip <- sign(x2-x1)
        x1 <- unit(x1,'native') ; x2 <- unit(x2,'native')
        xa <- x2 - unit(0.5*arrow.thick*flip,'char')
        xa <- if (flip > 0) max(xa,x1) else min(xa,x1)
        polygonGrob(unit.c(x1,xa,x2,xa,x1),
                    unit(margin+arrow.thick*c(0,0,0.5,1,1)+ofs,'char'),
                    gp=gpar(col=col, fill=col))
    }
    if (strand == '+') {
        a <- arrow(txStart+1, txEnd)
    } else {
        a <- arrow(txEnd, txStart+1)
    }
    xpos <- max(unit(min(txStart+1,txEnd),'native'), unit(1,'points'))
    xpos <- min(xpos, unit(1,'npc')-unit(1,'strwidth',name)-unit(1,'points'))
    ypos <- unit(margin+arrow.thick+text.sep+ofs,'char')
    s <- textGrob(name, xpos, ypos, c('left','bottom'))
    invisible(grobTree(a, s))
}

gene.track <-
    function(data, col='darkgreen', arrow.thick=0.4, margin=0.2, text.sep=0.1)
{
    if (identical(current.viewport()$xscale,c(0,1)))
        stop('track coordinates not initialized.  Use setup.tracks()')
    data$x.start <- convertX(unit(data$txStart,'native'),'npc',TRUE)
    data$x.end <- convertX(unit(data$txEnd,'native'),'npc',TRUE)
    data <- subset(data, (x.end > 0) & (x.start < 1))
    data <- data[order(data$txStart),]
    stru <- lapply(paste(data$geneSymbol,' ',sep=''),
                   unit, x=1, unit='strwidth')
    strw <- convertWidth(do.call('unit.c',stru), 'npc', TRUE)
    rmax <- c()
    xmin <- pmin(data$x.start,1-strw)
    xmax <- pmax(data$x.end,data$x.start+strw)
    for (i in 1:nrow(data)) {
        if (all(xmin[i] < rmax))
            rmax <- c(rmax,-Inf)
        rmax[which.max(xmin[i]-rmax)] <- xmax[i]
    }

    rmax <- rep(-Inf,length(rmax))
    names <- paste('row',1:length(rmax), sep='')
    gt <- grobTree(clipGrob())
    for (n in names) gt <- addGrob(gt, gTree(name=n))

    rh <- arrow.thick+margin+text.sep+1
    for (i in 1:nrow(data)) {
        rn <- which.max(xmin[i]-rmax)
        rmax[rn] <- xmax[i]
        g <- geneGrob(data$geneSymbol[i], data$strand[i],
                      data$txStart[i], data$txEnd[i], col,
                      rh*(rn-1), arrow.thick, margin, text.sep)
        gt <- addGrob(gt, g, gPath(names[rn]))
    }
    height <- unit(length(rmax)*rh+margin,'char')
    invisible(list(grob=gt, vp=list(height=height)))
}

manhattan.track <-
    function(x, y, threshold=-log10(5e-8), cex=0.25,
             col=c('#d0d0d0','#ff0000'), headroom=0,
             height=0.5*current.viewport()$width)
{
    if (identical(current.viewport()$xscale,c(0,1)))
        stop('track coordinates not initialized.  Use setup.tracks()')
    tck <- function(x) { p <- pretty(x) ; p[p >= x[1] & p <= x[2]] }
    w <- !is.na(x) & !is.na(y)
    x <- x[w] ; y <- y[w]
    ylim <- range(c(y,threshold))
    ylim <- ylim + c(-0.10,0.10)*diff(ylim)
    r <- 1 + headroom/(convertHeight(height,'lines',TRUE)-headroom)
    yscale <- ylim[1] + c(0,r*diff(ylim))

    ylab <- textGrob(expression(-log[10](pvalue)), x=unit(-2,'lines') -
                     unit(1,'strwidth',data=format(tck(ylim)))[1], rot=90)
    gt <- grobTree(yaxisGrob(tck(ylim), format(tck(ylim))),
                   yaxisGrob(tck(ylim), FALSE, FALSE), ylab,
                   linesGrob(y=rep(unit(threshold,'native'),2),
                             gp=gpar(col=col[1])),
                   pointsGrob(x, y, size=unit(cex,'char'),
                              gp=gpar(col=col[1+(y>threshold)])))
    invisible(list(grob=gt, vp=list(yscale=yscale, height=height)))
}
